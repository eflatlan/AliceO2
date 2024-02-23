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

/// \file   DigitsToRootSpec.cxx
/// \author Antonio Franco - INFN Bari
/// \version 1.0
/// \date 20 nov 2021
/// \brief Implementation of a data processor to produce Root File from Digits stream
///

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <array>
#include <functional>
#include <vector>
#include <algorithm>

#include "DPLUtils/DPLRawParser.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"

#include "TTree.h"
#include "TFile.h"

#include <gsl/span>

#include "Framework/DataRefUtils.h"
#include "Framework/InputSpec.h"
#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/Logger.h"
#include "Framework/InputRecordWalker.h"

#include "Headers/RAWDataHeader.h"
#include "DetectorsRaw/RDHUtils.h"

#include "HMPIDBase/Geo.h"
#include "HMPIDWorkflow/DigitsWriterSpec.h"

namespace o2
{
namespace hmpid
{

using namespace o2;
using namespace o2::framework;
using RDH = o2::header::RDHAny;

//=======================
void DigitsToRootTask::init(framework::InitContext& ic)
{
  LOG(info) << "[HMPID Write Root File From Digits stream - init()]";

  // get line command options
  mOutRootFileName = ic.options().get<std::string>("out-file");

  mTriggers.clear();
  mDigits.clear();

  TString filename = TString::Format("%s", mOutRootFileName.c_str());
  TString tit = TString::Format("HMPID Digits File Decoding");

  LOG(info) << "Create the ROOT file " << filename.Data();
  mfileOut.reset(new TFile(TString::Format("%s", filename.Data()), "RECREATE"));


// BranchDefinition<o2::dataformats::MCTruthContainer<o2::emcal::MCLabel>>{InputSpec{"emcaldigitlabels", "EMC", "DIGITSMCTR"}, "EMCALDigitMCTruth", mctruth ? 1 : 0})();

  mDigitTree.reset(new TTree("o2hmp", tit));
  mDigitTree->Branch("InteractionRecords", &mTriggers);
  mDigitTree->Branch("HMPIDDigits", &mDigits);
  
  if(mUseMC) {
    mDigitTree->Branch(mDigitMCTruthBranchName.c_str(), &mDigitLabels);
  }

  mExTimer.start();
  return;
}

void DigitsToRootTask::run(framework::ProcessingContext& pc)
{
  std::vector<o2::hmpid::Trigger> triggers;
  std::vector<o2::hmpid::Digit> digits;
  // std::vector<o2::hmpid::Digit> digits;
  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> digitLabels;
  
  
  for (auto const& ref : InputRecordWalker(pc.inputs())) {

    LOGP(info, "HMP CLASS DigitsToRootTask ");
    /*
     And what about MCLabels?
    */
    if (DataRefUtils::match(ref, {"check", ConcreteDataTypeMatcher{header::gDataOriginHMP, "INTRECORDS"}})) {
      triggers = pc.inputs().get<std::vector<o2::hmpid::Trigger>>(ref);
      LOG(info) << "We receive triggers =" << triggers.size();
    }
    if (DataRefUtils::match(ref, {"check", ConcreteDataTypeMatcher{header::gDataOriginHMP, "DIGITS"}})) {
      digits = pc.inputs().get<std::vector<o2::hmpid::Digit>>(ref);
      LOG(info) << "The size of the cluster-vector =" << digits.size();
    }
    
    
    
    //ef : we need to define the digitsMC-truth here ? 
    if (mUseMC) {
      if (DataRefUtils::match(ref, {"check", ConcreteDataTypeMatcher{header::gDataOriginHMP, "DIGITSMCTR"}})) {
        
        //std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> digitLabels
        // o2::dataformats::MCTruthContainer<o2::MCCompLabel>


        // we do this in the STEER? 
        digitLabels = pc.inputs().get<std::vector<o2::hmpid::Digit>>(ref);
        LOG(info) << "The size of the vector =" << digitLabels.size();
        LOGP(info, "HMP CLASS DigitsToRootTask got DIGITSMCTR");
      }
    }

    for (int i = 0; i < triggers.size(); i++) {
      LOG(info) << "Trigger Event     Orbit = " << triggers[i].getOrbit() << "  BC = " << triggers[i].getBc();
      int startDigitsIndex = mDigits.size();
      int numberOfDigits = 0;
      for (int j = triggers[i].getFirstEntry(); j <= triggers[i].getLastEntry(); j++) {
        mDigits.push_back(digits[j]); // append the cluster
        numberOfDigits++;


        // we do this in the STEER? 
        // if (mUseMC)
          // mDigitLabels.push_back(digitLabels[j]);
          // mDigitLabels.mergeAtBack(digitLabels[j]);?
      }

      mTriggers.push_back(triggers[i]);
      mTriggers.back().setDataRange(startDigitsIndex, numberOfDigits);
      // mDigitLabels.mergeAtBack(digitLabels);?

    }
  }

  mExTimer.stop();
  return;
}

void DigitsToRootTask::endOfStream(framework::EndOfStreamContext& ec)
{
  mExTimer.logMes("Received an End Of Stream !");
  LOG(info) << "The size of digits vector =" << mDigits.size();
  mDigitTree->Fill();
  mDigitTree->Write();
  mfileOut->Close();
  mExTimer.logMes("Register Tree ! ");
  return;
}

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getDigitsToRootSpec(std::string inputSpec, bool useMC)
{
  std::vector<o2::framework::InputSpec> inputs;


  // ???? ef :: TODO: why is ther clusters here?
  // ef; changed to digits
  inputs.emplace_back("digits", o2::header::gDataOriginHMP, "DIGITS", 0, Lifetime::Timeframe);
  inputs.emplace_back("intrecord", o2::header::gDataOriginHMP, "INTRECORDS", 0, Lifetime::Timeframe);


  //auto useMC = false; // ef not yet
  if(useMC) {
    inputs.emplace_back("hmpiddigitlabels", o2::header::gDataOriginHMP, "DIGITSMCTR", 0, Lifetime::Timeframe);
  } // ef: do as from the steer..


  
  std::vector<o2::framework::OutputSpec> outputs;

  return DataProcessorSpec{
    "HMP-DigitsToRoot",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<DigitsToRootTask>(useMC)},
    Options{{"out-file", VariantType::String, "hmpDigits.root", {"name of the output file"}}}};
}

} // namespace hmpid
} // end namespace o2
