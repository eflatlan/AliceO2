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
  mRec.reset(new o2::hmpid::Clusterer(mUseMC)); // ef: changed to smart-pointer
  mExTimer.start();
}

void DigitsToClustersTask::run(framework::ProcessingContext& pc)
{

  // bool mUseMC = true; // ef do inout tu fcn
  // outputs
  std::vector<o2::hmpid::Cluster> clusters;
  std::vector<o2::hmpid::Trigger> clusterTriggers;

  LOG(info) << "[HMPID DClusterization - run() ] Enter ...";

  auto triggers = pc.inputs().get<gsl::span<o2::hmpid::Trigger>>("intrecord");
  auto digits = pc.inputs().get<gsl::span<o2::hmpid::Digit>>("digits");

  // ef > added
  auto labelVector = std::make_shared<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>();
  if (mUseMC) {

    auto digitlabels = pc.inputs().get<o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>("hmpiddigitlabels");

    if (digitlabels == nullptr) {
      LOGP(error, "digitlabels nullptr");
      throw std::runtime_error("digitlabels was nullptr");
    }

    LOGP(info, "triggers {} : digits : {}", triggers.size(), digits.size());

    LOGP(info, "digitlabels of objs in truthArray : {}", digitlabels->getNElements());
    LOGP(info, "digitlabels of objs in headArray : {}", digitlabels->getIndexedSize());

    if (labelVector != nullptr) {
      *labelVector.get() = std::move(*digitlabels);
    } else {
      LOGP(error, "labelVector nullptr");
      throw std::runtime_error("labelVector was nullptr");
    }

    if (mClsLabels != nullptr) {
      mRec->setMCTruthContainer(mClsLabels.get());
    } else {
      LOGP(error, "mClsLabels nullptr");
      throw std::runtime_error("mClsLabels was nullptr");
    }
  }

  int i = 0;
  for (const auto& trig : triggers) {
    if (trig.getNumberOfObjects()) {

      LOGP(info, "[HMPID DClusterization - run() ]  trigger number {} number of digits : {} ", i, trig.getNumberOfObjects());

      gsl::span<const o2::hmpid::Digit> trigDigits{
        digits.data() + trig.getFirstEntry(),
        size_t(trig.getNumberOfObjects())};
      size_t clStart = clusters.size();

      if (mUseMC && labelVector != nullptr) {

        LOGP(info, "number of objs in truthArray : {}", labelVector->getNElements());
        LOGP(info, "number of objs in headArray : {}", labelVector->getIndexedSize());

        // clusLabels is set to Clusterer (mRec)
        // by : mRec->setMCTruthContainer(mClsLabels.get()

        mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, labelVector.get(), true);
        // mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, &(labelVector->at(i)), true);

      } else {
        // pass nullptr for mClustersArray if useMC == false
        mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, nullptr, true);
      }

      // ef > remove
      {
        auto timeA =
          o2::InteractionRecord::bc2ns(trig.getBc(), trig.getOrbit());
        int cnt = 0;
        int firstentry = trig.getFirstEntry();
        int lastEntry = trig.getLastEntry();
        LOGP(info,
             "START :DIGIT TRIGGER  entries {} first {}  lasrt {}  time {} ",
             trig.getNumberOfObjects(), firstentry, lastEntry, timeA / 1000.0f);

        LOGP(info, " bc {} orbit {} ", trig.getBc(), trig.getOrbit());
        LOGP(info, " entries {}", trig.getNumberOfObjects());
      }

      clusterTriggers.emplace_back(trig.getIr(), clStart, clusters.size() - clStart);

      // ef > remove
      auto t = clusterTriggers.back();
      {
        auto timeA = o2::InteractionRecord::bc2ns(t.getBc(), t.getOrbit());
        int cnt = 0;
        int firstentry = t.getFirstEntry();
        int lastEntry = t.getLastEntry();
        LOGP(info,
             "START : CLUSTER TRIGGER entries {} first {}  lasrt {}  time {} ",
             t.getNumberOfObjects(), firstentry, lastEntry, timeA / 1000.0f);

        LOGP(info, " bc {} orbit {} ", trig.getBc(), trig.getOrbit());
        LOGP(info, "entries {}", trig.getNumberOfObjects());
      }

    } else {
      LOGP(warn, "[HMPID DClusterization - run() ] trigger number {} had no digits", i);
    }
    i++;
  }

  LOGP(info, "Received {} triggers with {} digits -> {} triggers with {} clusters",
       triggers.size(), digits.size(), clusterTriggers.size(), clusters.size());
  mDigitsReceived += digits.size();
  mClustersReceived += clusters.size();

  pc.outputs().snapshot(o2::framework::Output{"HMP", "CLUSTERS", 0}, clusters);
  pc.outputs().snapshot(o2::framework::Output{"HMP", "INTRECORDS1", 0}, clusterTriggers);

  mExTimer.elapseMes("Clusterization of Digits received = " + std::to_string(mDigitsReceived));
  mExTimer.elapseMes("Clusterization of Clusters received = " + std::to_string(mClustersReceived));

  // ef: added mClsLabels
  if (mUseMC) {
    if (mClsLabels) {
      pc.outputs().snapshot(o2::framework::Output{"HMP", "CLUSTERSMCTR", 0}, *mClsLabels);
      LOGP(info, "[HMPID DigitsToCluster -  mcLabels size : headerArray {}; truthArray {}", mClsLabels->getIndexedSize(), mClsLabels->getNElements());
    }
  }
}

void DigitsToClustersTask::endOfStream(framework::EndOfStreamContext& ec)
{

  mExTimer.stop();
  mExTimer.logMes("End Clusterization !  digits = " +
                  std::to_string(mDigitsReceived));
}

//_______________________________________________________________________________________________
o2::framework::DataProcessorSpec
  getDigitsToClustersSpec(bool useMC)
{

  std::vector<o2::framework::InputSpec> inputs;

  inputs.emplace_back("digits", o2::header::gDataOriginHMP, "DIGITS", 0,
                      o2::framework::Lifetime::Timeframe);
  inputs.emplace_back("intrecord", o2::header::gDataOriginHMP, "INTRECORDS", 0,
                      o2::framework::Lifetime::Timeframe);

  // ef: added
  if (useMC) {
    inputs.emplace_back("hmpiddigitlabels", o2::header::gDataOriginHMP, "DIGITSMCTR", 0, Lifetime::Timeframe); // DIGITLBL == > DIGITSMCTR?
  }

  // define outputs
  std::vector<o2::framework::OutputSpec> outputs;

  // ef: added
  if (useMC) {
    outputs.emplace_back("HMP", "CLUSTERSMCTR", 0,
                         o2::framework::Lifetime::Timeframe);
  }

  outputs.emplace_back("HMP", "CLUSTERS", 0,
                       o2::framework::Lifetime::Timeframe);
  outputs.emplace_back("HMP", "INTRECORDS1", 0,
                       o2::framework::Lifetime::Timeframe);

  return DataProcessorSpec{
    "HMP-Clusterization", inputs, outputs,
    AlgorithmSpec{adaptFromTask<DigitsToClustersTask>(useMC)},
    Options{{"sigma-cut",
             VariantType::String,
             "",
             {"sigmas as comma separated list"}}}};
}

} // namespace hmpid
} // end namespace o2
