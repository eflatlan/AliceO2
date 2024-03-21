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
  /*clusters.clear();
  clusterTriggers.clear();*/

  // ef: added

  /*if (mUseMC) {
    mClsLabels.reset(new o2::dataformats::MCTruthContainer<o2::MCCompLabel>);
  }*/

  // bool mUseMC = false;// not yet

  /*
  auto labelVector = std::make_shared<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>>();
  if (mUseMC) {

    // ef FIX!
    //mClsLabels.clear();

    auto digitlabels = pc.inputs().get<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>*>("hmpiddigitlabels");
    *labelVector.get() = std::move(*digitlabels);
    //mRec->setMCTruthContainer(mClsLabels); ef do laters
  }*/

  // filled in HMPIDDigitizerSpec :
  // o2::dataformats::MCTruthContainer<o2::MCCompLabel> mLabels; // labels which get filled
  // std::vector<o2::hmpid::Trigger> mIntRecord;

  auto triggers = pc.inputs().get<gsl::span<o2::hmpid::Trigger>>("intrecord");
  auto digits = pc.inputs().get<gsl::span<o2::hmpid::Digit>>("digits");

  // auto labelVector = std::make_shared<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>>();

  auto labelVector = std::make_shared<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>();
  if (mUseMC) {

    LOGP(info, "Trying to acces digitLabels");

    // ef: I dont understand which to use :
    // auto digitlabels = pc.inputs().get<o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>("hmpiddigitlabels");

    // ef: if based on other.
    auto digitlabels = pc.inputs().get<o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>("hmpiddigitlabels");

    if (digitlabels == nullptr) {
      LOGP(info, "digitlabels nullptr");
    }

    // B) ef: if based on TOF:
    // auto digitlabels = pc.inputs().get<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>*>("hmpiddigitlabels");
    // LOGP(info, "digitlabels size : {}", digitlabels->getSize());

    LOGP(info, "triggers {} : digits : {}", triggers.size(), digits.size());

    LOGP(info, "digitlabels of objs in truthArray : {}", digitlabels->getNElements());
    LOGP(info, "digitlabels of objs in headArray : {}", digitlabels->getIndexedSize());

    if (digitlabels != nullptr)
      *labelVector.get() = std::move(*digitlabels);

    if (mClsLabels != nullptr) {
      mRec->setMCTruthContainer(mClsLabels.get());
      LOGP(info, "mClsLabels was not nullpttr");
    } else {
      LOGP(info, "mClsLabels was nullptr");
    }
  }

  int i = 0;
  LOGP(info, "loop of triggers");
  for (const auto& trig : triggers) {
    if (trig.getNumberOfObjects()) {

      LOGP(info, "[HMPID DClusterization - run() ] trig.numObjects : {} <<< calling clusterer->Dig2Clu", trig.getNumberOfObjects());

      gsl::span<const o2::hmpid::Digit> trigDigits{
        digits.data() + trig.getFirstEntry(),
        size_t(trig.getNumberOfObjects())};
      size_t clStart = clusters.size();
      LOGP(info, "[HMPID DClusterization  clStart {}", clStart);

      // ef; taken from TOF, does it work?
      // &(labelVector->at(i)) ?

      // A)  labelVector->at(i) : if we use vector<MCTruthContainer>
      // B) labelVector : if we use MCTruthContainer

      // as of now we store triggers in vector, and
      // mcTruth for digits as MCTruthContainer<o2::MCCompLabel> :
      // https://github.com/AliceO2Group/AliceO2/blob/25fed0222034939656422b6c0d727dd9cad1983b/Steer/DigitizerWorkflow/src/HMPIDDigitizerSpec.cxx#L92

      //

      // ef : I am not sure about labelVector, should we not do as labelVector->at(i) ?
      // 		to sort this based on each trigger?
      //    I anyway changed to store the label of the digits in Digitizer
      //    So it should be able to seperate based on tjis
      //    (Also seperate which event the digit corresponds to )

      if (mUseMC && labelVector != nullptr) {
        LOGP(info, "[HMPID DClusterization mUseMC {}", mUseMC);

        // TOF mRec.process(mReader, mClustersArray, &(labelVector->at(i)));
        /*if(labelVector==nullptr) {
          LOGP(info, "labelVector was nullptr");
        } else {
          LOGP(info, "labelVector of size {}", labelVector->size());
        }*/

        /*if(digitlabels==nullptr) {
        }*/

        LOGP(info, "IP args to Dig2Clu : ");
        LOGP(info, " clusters Size {}", clusters.size());
        LOGP(info, " trigDigits Size {}", trigDigits.size());

        // add protection against this ?
        /*if(i > labelVector->size()) {
          LOGP(info, "labelVector->size() {}", labelVector->size());
        }*/

        LOGP(info, "number of objs in truthArray : {}", labelVector->getNElements());
        LOGP(info, "number of objs in headArray : {}", labelVector->getIndexedSize());

        // auto labelObj = labelVector->getLabels(i);

        // clusLabels is set to Clusterer (mRec)
        // by : mRec->setMCTruthContainer(mClsLabels.get()

        // ef :should this not be vectors?
        mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, labelVector.get(), true);

        // mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, &(labelVector->at(i)), true);

        LOGP(info, "[HMPID DClusterization exit :: mRec->Dig2Clu");

      } else {
        // mRec.process(mReader, mClustersArray, nullptr);
        mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, nullptr, true);
      }

      LOG(info) << "[HMPID DClusterization - return from dig2clu";

      // if(clusters.back().dig(0) == nullptr) {Printf("DigtisToClusterSpec:: dig was nullptr!!");}
      clusterTriggers.emplace_back(trig.getIr(), clStart, clusters.size() - clStart);
    } else {

      LOGP(info, "[HMPID DClusterization - run() ] trig.numObjects 0");
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
      LOGP(info, "[HMPID DigitsToCluster -  mcLabels size : headerArray {}; truthArray {}", mClsLabels->getIndexedSize(),  mClsLabels->getNElements()); 
    }
  }



  /* ef: should it be cleared?
  if(mClsLabels != nullptr) {
    LOGP(info, "mClsLabels->clear()");
    mClsLabels->clear();
  } else {
    LOGP(info, "mClsLabels was nullptr");
  }*/
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
    inputs.emplace_back("hmpiddigitlabels", o2::header::gDataOriginHMP, "DIGITSMCTR", 0, Lifetime::Timeframe);
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
