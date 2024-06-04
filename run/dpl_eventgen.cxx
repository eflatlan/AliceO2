// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "../Framework/Core/src/ArrowSupport.h"
#include "Framework/AnalysisTask.h"
#include "Monitoring/Monitoring.h"
#include "Framework/CommonDataProcessors.h"
#include "SimulationDataFormat/MCTrack.h"
#include "Framework/runDataProcessing.h"
#include <Generators/GeneratorService.h>
#include <CommonUtils/ConfigurableParam.h>
#include <CommonUtils/RngHelper.h>
#include <TStopwatch.h> // simple timer from ROOT
#include <DetectorsCommonDataFormats/DetectorNameConf.h>
#include <DetectorsBase/Detector.h>
#include <TFile.h>
#include <memory>

using namespace o2::framework;

struct GeneratorTask {
  // For readability to indicate where counting certain things (such as events or timeframes) should be of the same order of magnitude
  typedef uint64_t GenCount;
  Configurable<std::string> generator{"generator", "boxgen", "Name of generator"};
  Configurable<GenCount> eventNum{"nEvents", 1, "Number of events"};
  Configurable<std::string> trigger{"trigger", "", "Trigger type"}; //
  Configurable<std::string> iniFile{"configFile", "", "INI file containing configurable parameters"};
  Configurable<std::string> params{"configKeyValues", "", "configurable params - configuring event generation internals"};
  Configurable<long> seed{"seed", 0, "(TRandom) Seed"};
  Configurable<int> aggregate{"aggregate-timeframe", 300, "Number of events to put in a timeframe"};
  Configurable<std::string> vtxModeArg{"vertexMode", "kDiamondParam", "Where the beam-spot vertex should come from. Must be one of kNoVertex, kDiamondParam, kCCDB"};
  Configurable<int64_t> ttl{"time-limit", -1, "Maximum run time limit in seconds (default no limit)"};
  Configurable<std::string> outputPrefix{"output", "", "Optional prefix for kinematics files written on disc. If non-empty, files <prefix>_Kine.root + <prefix>_MCHeader.root will be created."};
  GenCount nEvents = 0;
  GenCount eventCounter = 0;
  GenCount tfCounter = 0;
  std::unique_ptr<TFile> outfile{};
  std::unique_ptr<TTree> outtree{};

  // a pointer because object should only be constructed in the device (not during DPL workflow setup)
  std::unique_ptr<o2::eventgen::GeneratorService> genservice;
  TStopwatch timer;

  void init(o2::framework::InitContext& /*ic*/)
  {
    genservice.reset(new o2::eventgen::GeneratorService);
    o2::utils::RngHelper::setGRandomSeed(seed);
    nEvents = eventNum;
    // helper to parse vertex option; returns true if parsing ok, false if failure
    o2::conf::VertexMode vtxmode;
    if (!(o2::conf::SimConfig::parseVertexModeString(vtxModeArg, vtxmode))) {
      LOG(error) << "Could not parse vtxMode";
    }

    // update config key params
    o2::conf::ConfigurableParam::updateFromFile(iniFile);
    o2::conf::ConfigurableParam::updateFromString((std::string)params);
    // initialize the service
    if (vtxmode == o2::conf::VertexMode::kDiamondParam) {
      genservice->initService(generator, trigger, o2::eventgen::DiamondParamVertexOption());
    } else if (vtxmode == o2::conf::VertexMode::kNoVertex) {
      genservice->initService(generator, trigger, o2::eventgen::NoVertexOption());
    } else if (vtxmode == o2::conf::VertexMode::kCCDB) {
      LOG(warn) << "Not yet supported. This needs definition of a timestamp and fetching of the MeanVertex CCDB object";
    }
    timer.Start();

    if (outputPrefix->size() > 0 && !outfile.get()) {
      auto kineoutfilename = o2::base::NameConf::getMCKinematicsFileName(outputPrefix->c_str());
      outfile.reset(new TFile(kineoutfilename.c_str(), "RECREATE"));
      outtree.reset(new TTree("o2sim", "o2sim"));
    }
  }

  void run(o2::framework::ProcessingContext& pc)
  {
    std::vector<o2::MCTrack> mctracks;
    o2::dataformats::MCEventHeader mcheader;
    auto mctrack_ptr = &mctracks;
    if (outfile.get()) {
      auto br = o2::base::getOrMakeBranch(*outtree, "MCTrack", &mctrack_ptr);
      br->SetAddress(&mctrack_ptr);
    }

    for (auto i = 0; i < std::min((GenCount)aggregate, nEvents - eventCounter); ++i) {
      mctracks.clear();
      genservice->generateEvent_MCTracks(mctracks, mcheader);
      pc.outputs().snapshot(Output{"MC", "MCHEADER", 0}, mcheader);
      pc.outputs().snapshot(Output{"MC", "MCTRACKS", 0}, mctracks);
      ++eventCounter;

      if (outfile.get() && outtree.get()) {
        outtree->Fill();
      }
    }

    // report number of TFs injected for the rate limiter to work
    ++tfCounter;
    pc.services().get<o2::monitoring::Monitoring>().send(o2::monitoring::Metric{(uint64_t)tfCounter, "df-sent"}.addTag(o2::monitoring::tags::Key::Subsystem, o2::monitoring::tags::Value::DPL));
    bool time_expired = false;
    if (ttl > 0) {
      timer.Stop();
      time_expired = timer.RealTime() > ttl;
      timer.Start(false);
      if (time_expired) {
        LOG(info) << "TTL expired after " << eventCounter << " events ... sending end-of-stream";
      }
    }
    if (eventCounter >= nEvents || time_expired) {
      pc.services().get<ControlService>().endOfStream();
      pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);

      // write out data to disc if asked
      if (outfile.get()) {
        outtree->SetEntries(eventCounter);
        outtree->Write();
        outfile->Close();
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto spec = adaptAnalysisTask<GeneratorTask>(cfgc);
  spec.outputs.emplace_back("MC", "MCHEADER", 0, Lifetime::Timeframe);
  spec.outputs.emplace_back("MC", "MCTRACKS", 0, Lifetime::Timeframe);
  spec.requiredServices.push_back(o2::framework::ArrowSupport::arrowBackendSpec());
  spec.algorithm = CommonDataProcessors::wrapWithRateLimiting(spec.algorithm);
  return {spec};
}
