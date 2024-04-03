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

/// \file   ClusterReaderSpec.cxx
/// \author Annalisa Mastroserio - INFN Bari
/// \version 1.0
/// \date 22 Jun 2022
/// \brief Implementation of a data processor to read Cluster tree and provide the array for further usage
///

#include "HMPIDWorkflow/ClustersReaderSpec.h"

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <array>
#include <functional>
#include <vector>

#include "CommonUtils/StringUtils.h" // o2::utils::Str

#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Logger.h"
#include "Framework/DataRefUtils.h"
#include "Framework/InputRecordWalker.h"

#include "Headers/RAWDataHeader.h"
#include "DetectorsRaw/RDHUtils.h"
#include "DPLUtils/DPLRawParser.h"

namespace o2
{
namespace hmpid
{

using namespace o2;
/// using namespace o2::header;
using namespace o2::framework;
using RDH = o2::header::RDHAny;

//
void ClusterReaderTask::init(framework::InitContext& ic)
{
  LOG(info) << "[HMPID Cluster reader - init() ] ";
  mClustersReceived = 0;
  // Build the file name
  const auto filename = o2::utils::Str::concat_string(
    o2::utils::Str::rectifyDirectory(
      ic.options().get<std::string>("input-dir")),
    ic.options().get<std::string>("hmpid-cluster-infile" /*"qc-hmpid-clusters"*/));
  initFileIn(filename);
}

void ClusterReaderTask::run(ProcessingContext& pc)
{

  auto ent = mTree->GetReadEntry() + 1;
  assert(ent < mTree->GetEntries()); // this should not happen
  mTree->GetEntry(ent);

  pc.outputs().snapshot(Output{"HMP", "CLUSTERS", 0}, mClustersFromFile);
  pc.outputs().snapshot(Output{"HMP", "INTRECORDS1", 0}, mClusterTriggersFromFile);
  mClustersReceived += mClustersFromFile.size();
  LOG(info) << "[HMPID ClusterReader - run() ] clusters  = " << mClustersFromFile.size();

  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
    pc.services().get<ControlService>().endOfStream();
    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    mExTimer.stop();
    mExTimer.logMes("End ClusterReader !  clusters = " +
                    std::to_string(mClustersReceived));
  }


    int tnum = 0;

    LOGP(info, "triggers > {}", mClusterTriggersFromFile.size());
    for (const auto trig : *mClusterTriggersFromFilePtr) {

      auto timeA = o2::InteractionRecord::bc2ns(trig.getBc(), trig.getOrbit());
      int cnt = 0;
      int firstentry = trig.getFirstEntry(); int lastEntry = trig.getLastEntry();
      LOGP(info, "START : trigger number {} : entries {} first {}  lasrt {}  time {} ",tnum, trig.getNumberOfObjects(),  firstentry, lastEntry, timeA / 1000.0f);

      LOGP(info, " bc {} orbit {} ", trig.getBc(), trig.getOrbit());
      LOGP(info, "end{} entries {}", tnum, trig.getNumberOfObjects());
      tnum++;
    }





}

void ClusterReaderTask::initFileIn(const std::string& filename)
{
  // Create the TFIle
  mTree.reset(nullptr);
  mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());

  if (!mFile->IsOpen() || mFile == nullptr) {
    LOG(error) << "HMPID ClusterReaderTask::init() : Did not find any Clusters file " << filename.c_str() << " file !";
    throw std::runtime_error("cannot open input clusters file");
  } else {
    LOG(info) << "HMPID ClusterReaderTask::init() : Found Clusters file " << filename.c_str();
  }

  if ((TTree*)mFile->Get("o2hmp")) {
    mTree.reset((TTree*)mFile->Get("o2hmp"));
  } else if ((TTree*)mFile->Get("o2sim")) {
    mTree.reset((TTree*)mFile->Get("o2sim"));
  } else {
    LOG(error)
      << "HMPID ClusterReaderTask::init() : Did not find either Tree o2sim or o2hmp tree in "
      << filename.c_str();
    throw std::runtime_error(
      "HMPID ClusterReaderTask::init() : Did not find "
      "o2sim file in clusters tree");
  }

  if (mTree->GetBranchStatus("HMPIDClusters") == 1) {
    mTree->SetBranchAddress("HMPIDClusters", &mClustersFromFilePtr);
  } else if (mTree->GetBranchStatus("HMPIDclusters") == 1) {
    mTree->SetBranchAddress("HMPIDclusters", &mClustersFromFilePtr);
  } else {
    LOG(error)
      << "HMPID ClusterReaderTask::init() : Did not find Branch in "
      << filename.c_str();
    throw std::runtime_error(
      "HMPID ClusterReaderTask::init() : Did not find Branch HMPIDClusters in clusters tree");
    }

    int tnum = 0;

    LOGP(info, "triggers > {}", mClusterTriggersFromFile.size());
    for (const auto trig : *mClusterTriggersFromFilePtr) {

      auto timeA = o2::InteractionRecord::bc2ns(trig.getBc(), trig.getOrbit());
      int cnt = 0;
      int firstentry = trig.getFirstEntry(); int lastEntry = trig.getLastEntry();
      LOGP(info, "START : trigger number {} : entries {} first {}  lasrt {}  time {} ",tnum, trig.getNumberOfObjects(),  firstentry, lastEntry, timeA / 1000.0f);

      LOGP(info, " bc {} orbit {} ", trig.getBc(), trig.getOrbit());


      /*std::vector<int> eventLabels;

      for (int i = trig.getFirstEntry(); i <= trig.getLastEntry(); i++) {

        /*****insert from digitsreder* / //

        if (i < mLabels.getIndexedSize() && i < mClustersFromFile.size()) {
          bool isLabelEventSame = true;
          const auto& labels = mLabels.getLabels(i);
          int prevEventLabel = 0;

          if (labels.size() > 0) {
            prevEventLabel = labels[0].getEventID();
            eventLabels.push_back(prevEventLabel);
          }
          int lblNum = 0;
          for (const auto& label : labels) {

            if (label.getEventID() != prevEventLabel) {
              eventLabels.push_back(label.getEventID());
              isLabelEventSame = false;
              // LOGP(info, "trigger number {} lblNum {} : event from labelEventId changed!", tnum, lblNum);
              // LOGP(info, "digit number {}, digEventNum {} labelEventId {} prevEventLabel {}", i, mClustersFromFile[i].getEventNumber(), label.getEventID(), prevEventLabel);
            }
            lblNum++;

            prevEventLabel = label.getEventID();
          }
        }
      } 
      */

 

      /*****insert from digitsreder*/ //

      LOGP(info, "end{} entries {}", tnum, trig.getNumberOfObjects());
      tnum++;
    }


  mTree->SetBranchAddress("InteractionRecords", &mClusterTriggersFromFilePtr);
  mTree->Print("toponly");
}

//_________________________________________________________________________________________________

o2::framework::DataProcessorSpec getClusterReaderSpec()
{

  std::vector<o2::framework::OutputSpec> outputs;
  outputs.emplace_back("HMP", "CLUSTERS", 0, o2::framework::Lifetime::Timeframe);
  outputs.emplace_back("HMP", "INTRECORDS1", 0, o2::framework::Lifetime::Timeframe);

  return DataProcessorSpec{
    "HMP-ClusterReader",
    Inputs{},
    outputs,
    AlgorithmSpec{adaptFromTask<ClusterReaderTask>()},
    Options{{"hmpid-cluster-infile" /*"qc-hmpid-clusters"*/, VariantType::String, "hmpidclusters.root", {"Name of the input file with clusters"}},
            {"input-dir", VariantType::String, "./", {"Input directory"}}}};
}

} // namespace hmpid
} // end namespace o2
