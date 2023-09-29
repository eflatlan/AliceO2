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

/// \file DigitsToClustersSpec.cxx
/// \brief Implementation of clusterization for HMPID; read upstream/from file write upstream/to file

#include "HMPIDWorkflow/DigitsToClustersSpec.h"
#include "HMPIDWorkflow/DigitsReaderSpec.h"
#include "HMPIDWorkflow/ClustersWriterSpec.h"

#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataRefUtils.h"
#include "Framework/InputRecordWalker.h"
#include "Framework/Lifetime.h"
#include "Framework/Logger.h"
#include "Framework/Output.h"

#include "DPLUtils/DPLRawParser.h"
#include "DetectorsRaw/RDHUtils.h"
#include "Headers/RAWDataHeader.h"

#include "CommonUtils/NameConf.h" // o2::utils::Str

namespace o2
{
namespace hmpid
{

using namespace o2;
using namespace o2::header;
using namespace o2::framework;
using RDH = o2::header::RDHAny;

// Splits a string in float array for string delimiter, TODO: Move this in a
// HMPID common library
void DigitsToClustersTask::strToFloatsSplit(std::string s,
                                            std::string delimiter, float* res,
                                            int maxElem)
{
  int index = 0;
  size_t pos_start = 0;
  size_t pos_end;
  size_t delim_len = delimiter.length();
  std::string token;
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res[index++] = std::stof(token);
    if (index == maxElem) {
      return;
    }
  }
  res[index++] = (std::stof(s.substr(pos_start)));
  return;
}

//=======================
//
void DigitsToClustersTask::init(framework::InitContext& ic)
{
  mSigmaCutPar = ic.options().get<std::string>("sigma-cut");

  if (mSigmaCutPar != "") {
    strToFloatsSplit(mSigmaCutPar, ",", mSigmaCut, 7);
  }

  mDigitsReceived = 0, mClustersReceived = 0;

  mRec.reset(new o2::hmpid::Clusterer()); // ef: changed to smart-pointer

  mExTimer.start();
}

void DigitsToClustersTask::run(framework::ProcessingContext& pc)
{

	bool mUseMC = false; // ef do inout tu fcn
  // outputs
  std::vector<o2::hmpid::Cluster> clusters;
  std::vector<o2::hmpid::Trigger> clusterTriggers;
  std::vector<o2::hmpid::Topology> topVectorVector;

  LOG(info) << "[HMPID DClusterization - run() ] Enter ...";
  /*clusters.clear();
  clusterTriggers.clear();*/ 

  // ef: added

  /*if (mUseMC) {
    mClsLabels.reset(new o2::dataformats::MCTruthContainer<o2::MCCompLabel>);
  }*/ 

  //bool mUseMC = false;// not yet

  /*
  auto labelvector = std::make_shared<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>>();
  if (mUseMC) {

    // ef FIX!
    //mClsLabels.clear();
    
    auto digitlabels = pc.inputs().get<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>*>("tofdigitlabels");
    *labelvector.get() = std::move(*digitlabels);
    //mRec->setMCTruthContainer(mClsLabels); ef do laters
  }*/


  auto triggers = pc.inputs().get<gsl::span<o2::hmpid::Trigger>>("intrecord");
  auto digits = pc.inputs().get<gsl::span<o2::hmpid::Digit>>("digits");


  int i = 0;
  for (const auto& trig : triggers) {
    if (trig.getNumberOfObjects()) {

  		LOG(info) << "[HMPID DClusterization - run() ] nObjects <<< calling clusterer->Dig2Clu" << trig.getNumberOfObjects();
      gsl::span<const o2::hmpid::Digit> trigDigits{
        digits.data() + trig.getFirstEntry(),
        size_t(trig.getNumberOfObjects())};
      size_t clStart = clusters.size();

	
      if (mUseMC) {
        // TOF mClusterer.process(mReader, mClustersArray, &(labelvector->at(i)));
        //mRec->Dig2Clu(trigDigits, clusters, topVectorVector, mSigmaCut, &(labelvector->at(i)),true); ef do this later
        mRec->Dig2Clu(trigDigits, clusters, topVectorVector, mSigmaCut, nullptr, true);
      } else {
        // mClusterer.process(mReader, mClustersArray, nullptr);
        mRec->Dig2Clu(trigDigits, clusters, topVectorVector, mSigmaCut, nullptr, true);
      }


       LOG(info) << "[HMPID DClusterization - return from dig2clu";

      //if(clusters.back().dig(0) == nullptr) {Printf("DigtisToClusterSpec:: dig was nullptr!!");}
      clusterTriggers.emplace_back(trig.getIr(), clStart, clusters.size() - clStart);
    }
    i++;
  }
  LOGP(info, "Received {} triggers with {} digits -> {} triggers with {} clusters",
       triggers.size(), digits.size(), clusterTriggers.size(), clusters.size());
  mDigitsReceived += digits.size();
  mClustersReceived += clusters.size();


 /*ef: FIX
 if (mUseMC) {
      pc.outputs().snapshot(o2::framework::Output{"HMP", "CLUSTERSMCTR", 0, o2::framework::Lifetime::Timeframe}, mClsLabels);
  } */

  pc.outputs().snapshot(o2::framework::Output{"HMP", "CLUSTERS", 0, o2::framework::Lifetime::Timeframe}, clusters);

  pc.outputs().snapshot(o2::framework::Output{"HMP", "DIGITTOPOLOGY", 0, o2::framework::Lifetime::Timeframe}, topVectorVector);

  pc.outputs().snapshot(o2::framework::Output{"HMP", "INTRECORDS1", 0, o2::framework::Lifetime::Timeframe}, clusterTriggers);

  mExTimer.elapseMes("Clusterization of Digits received = " + std::to_string(mDigitsReceived));
  mExTimer.elapseMes("Clusterization of Clusters received = " + std::to_string(mClustersReceived));
}

void DigitsToClustersTask::endOfStream(framework::EndOfStreamContext& ec)
{

  mExTimer.stop();
  mExTimer.logMes("End Clusterization !  digits = " +
                  std::to_string(mDigitsReceived));
}

//_______________________________________________________________________________________________
o2::framework::DataProcessorSpec
  getDigitsToClustersSpec()

{

  std::vector<o2::framework::InputSpec> inputs;

  inputs.emplace_back("digits", o2::header::gDataOriginHMP, "DIGITS", 0,
                      o2::framework::Lifetime::Timeframe);
  inputs.emplace_back("intrecord", o2::header::gDataOriginHMP, "INTRECORDS", 0,
                      o2::framework::Lifetime::Timeframe);

  bool mUseMC = false; // ef do later

  if (mUseMC) {
    inputs.emplace_back("tofdigitlabels", o2::header::gDataOriginTOF, "DIGITSMCTR", 0, Lifetime::Timeframe);
  }

  // define outputs
  std::vector<o2::framework::OutputSpec> outputs;


 // ef: FIX
  if (mUseMC) {
  outputs.emplace_back("HMP", "CLUSTERSMCTR", 0,
                       o2::framework::Lifetime::Timeframe);
  }

  outputs.emplace_back("HMP", "CLUSTERS", 0,
                       o2::framework::Lifetime::Timeframe);
  outputs.emplace_back("HMP", "DIGITTOPOLOGY", 0,
                       o2::framework::Lifetime::Timeframe);
  outputs.emplace_back("HMP", "INTRECORDS1", 0,
                       o2::framework::Lifetime::Timeframe);

  return DataProcessorSpec{
    "HMP-Clusterization", inputs, outputs,
    AlgorithmSpec{adaptFromTask<DigitsToClustersTask>()},
    Options{{"sigma-cut",
             VariantType::String,
             "",
             {"sigmas as comma separated list"}}}};
}

} // namespace hmpid
} // end namespace o2
