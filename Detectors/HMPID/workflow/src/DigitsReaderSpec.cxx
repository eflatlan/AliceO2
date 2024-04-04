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

/// @file   DigitReaderSpec.cxx

#include <vector>

#include <TTree.h>

#include "Framework/ConfigParamRegistry.h"

#include "Framework/ControlService.h"

#include "Framework/Logger.h"

#include "SimulationDataFormat/MCTruthContainer.h"

#include "CommonUtils/NameConf.h"

#include "HMPIDWorkflow/DigitsReaderSpec.h"

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

#include "Framework/WorkflowSpec.h"

#include "Framework/Logger.h"

#include "Framework/DataRefUtils.h"

#include "Framework/InputRecordWalker.h"

#include "Headers/RAWDataHeader.h"

#include "DetectorsRaw/RDHUtils.h"

#include "DPLUtils/DPLRawParser.h"

#include "HMPIDBase/Geo.h"

using namespace o2::framework;

using namespace o2;

using namespace o2::header;

namespace o2

{

namespace hmpid

{

void DigitReader::init(InitContext& ic)

{

  auto filename = o2::utils::Str::concat_string(o2::utils::Str::rectifyDirectory(ic.options().get<std::string>("input-dir")),

                                                ic.options().get<std::string>("hmpid-digit-infile"));

  mFile.reset(TFile::Open(filename.c_str()));

  LOG(info) << "HMPID DigitWriterSpec::init() : Trying to read File : " << filename.c_str();

  mDigitsReceived = 0;

  if (!mFile->IsOpen()) {

    LOG(error) << "HMPID DigitWriterSpec::init() : Did not find any digits file " << filename.c_str() << " file !";

    throw std::runtime_error("cannot open input digits file");
  }

  if ((TTree*)mFile->Get("o2sim") != nullptr) {

    mTreeDig.reset((TTree*)mFile->Get("o2sim"));

  } else if ((TTree*)mFile->Get("o2hmp") != nullptr) {

    mTreeDig.reset((TTree*)mFile->Get("o2hmp"));

  } else {

    LOG(error) << "Did not find o2hmp tree in " << filename.c_str();

    throw std::runtime_error("Did Not find Any correct Tree in HMPID Digits File");
  }

  if (!mTreeDig) {

    LOG(error) << "Did not find o2hmp tree in " << filename.c_str();

    throw std::runtime_error("Did Not find Any correct Tree in HMPID Digits File");

  } /*else {

    LOG(info) << "HMPID DigitWriterSpec::init() : Reading From Branch  o2hmp" << File : " << filename.c_str()";

  } */
}

void DigitReader::run(ProcessingContext& pc)

{

  std::vector<o2::hmpid::Digit> mDigitsFromFile, *mDigitsFromFilePtr = &mDigitsFromFile;

  std::vector<o2::hmpid::Trigger> mTriggersFromFile, *mTriggersFromFilePtr = &mTriggersFromFile;

  /*  */

  if (mTreeDig->GetBranchStatus("HMPDigit")) {

    mTreeDig->SetBranchAddress("HMPDigit", &mDigitsFromFilePtr);

  } else if (mTreeDig->GetBranchStatus("HMPIDDigits")) {

    mTreeDig->SetBranchAddress("HMPIDDigits", &mDigitsFromFilePtr);

  } else {

    LOG(error)

      << "HMPID DigitWriterSpec::init() : Did not find any branch for Digits";

    throw std::runtime_error("Did Not find Any correct Branch for Digits in HMPID Digits File");
  }

  if (mTreeDig->GetBranchStatus("InteractionRecords")) {

    mTreeDig->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);

  } else {

    LOG(error)

      << "HMPID DigitWriterSpec::init() : Did not find  branch for Triggers";

    throw std::runtime_error("Did Not find Branch For triggers in HMPID Digits File");
  }

  mTreeDig->Print("toponly");

  if (mTreeDig->GetBranchStatus("HMPIDDigitMCTruth")) {

    mTreeDig->SetBranchAddress("HMPIDDigitMCTruth", &mPlabels);
  }

  auto ent = mTreeDig->GetReadEntry() + 1;

  assert(ent < mTreeDig->GetEntries()); // this should not happen

  mTreeDig->GetEntry(ent);

  pc.outputs().snapshot(Output{"HMP", "DIGITS", 0}, mDigitsFromFile);

  pc.outputs().snapshot(Output{"HMP", "INTRECORDS", 0}, mTriggersFromFile);

  if (mUseMC) {

      pc.outputs().snapshot(Output{"HMP", "DIGITSMCTR", 0}, mLabels);
  }

  // ef remove latere :
  if (mVerbose) {

      int tnum = 0;

      for (const auto trig : mTriggersFromFile) {

        auto timeA =
            o2::InteractionRecord::bc2ns(trig.getBc(), trig.getOrbit());
        int firstentry = trig.getFirstEntry();
        int lastEntry = trig.getLastEntry();
        LOGP(info,
             "START : trigger number {} : entries {} first {}  lasrt {}  time "
             "{} ",
             tnum, trig.getNumberOfObjects(), firstentry, lastEntry,
             timeA / 1000.0f);
        //  LOGP(info, "START : trigger number {} : entries {} first {}  lasrt
        //  {}  time {} ",tnum, trig.getNumberOfObjects(),  firstentry,
        //  lastEntry, timeA / 1000.0f);

        LOGP(info, " bc {} orbit {} ", trig.getBc(), trig.getOrbit());
        int cnt = 0;
				if (mUseMC) {
		      std::vector<int> eventLabels;

		      for (int i = trig.getFirstEntry(); i <= trig.getLastEntry(); i++) {

		        if (i < mLabels.getIndexedSize() && i < mDigitsFromFile.size()) {

		          bool isLabelEventSame = true;
		          const auto& labels = mLabels.getLabels(i);
		          int prevEventLabel;
		          if(labels.size() > 0) {
		            prevEventLabel = labels[0].getEventID();
		            eventLabels.push_back(prevEventLabel);
		          }
		          int lblNum = 0;
		          for (const auto& label : labels) {

		            if(label.getEventID() != prevEventLabel) {

		              eventLabels.push_back(label.getEventID());
		              isLabelEventSame = false;

		              // LOGP(info, "trigger number {} lblNum {} : event from labelEventId changed!", tnum, lblNum);

		              // LOGP(info, "digit number {}, digEventNum {} labelEventId {} prevEventLabel {}", i, mDigitsFromFile[i].getEventNumber(), label.getEventID(), prevEventLabel);
		            }

		            lblNum++;

		            prevEventLabel = label.getEventID();

		          }
		        }
		        else {

		          LOGP(info, "out of range {} > numLabels {}", i, mLabels.getIndexedSize());
		        }
		      }

		      LOGP(info, "trigger number {} : entries {}", tnum,
		           trig.getNumberOfObjects());

		      LOGP(info, "\ndifferent labels from eventLabels {} :::", eventLabels.size());
		      std::vector<int> sortedVec = eventLabels;

		      std::sort(sortedVec.begin(), sortedVec.end());

		      std::cout << "eventLabels values: ";
		      for (size_t i = 0; i < sortedVec.size(); ++i) {
		        if (i == sortedVec.size() - 1 || sortedVec[i] != sortedVec[i + 1]) {
		          std::cout << sortedVec[i] << " , ";
		        }
		      }
				}
        LOGP(info, "cnt {} entries {}", cnt, trig.getNumberOfObjects());

        tnum++;
				
      }
  }

  mDigitsReceived += mDigitsFromFile.size();

  LOG(info) << "[HMPID DigitsReader - run() ] digits  = " << mDigitsFromFile.size();

  if (mTreeDig->GetReadEntry() + 1 >= mTreeDig->GetEntries()) {

    pc.services().get<ControlService>().endOfStream();

    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);

    mExTimer.stop();

    if(mUseMC) {

      LOGP(info, "[HMPID DigitsReader - with useMC : mcLabels size : headerArray {}; truthArray {}", mLabels.getIndexedSize(),  mLabels.getNElements());
    }

    mExTimer.logMes("End DigitsReader !  digits = " +

                    std::to_string(mDigitsReceived));
  }
}

DataProcessorSpec getDigitsReaderSpec(bool useMC, bool verbose)
{

  std::vector<OutputSpec> outputs;

  outputs.emplace_back("HMP", "DIGITS", 0, o2::framework::Lifetime::Timeframe);

  outputs.emplace_back("HMP", "INTRECORDS", 0, o2::framework::Lifetime::Timeframe);

  if (useMC) {
    outputs.emplace_back("HMP", "DIGITSMCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{

    "HMP-DigitReader",

    Inputs{},

    outputs,

    AlgorithmSpec{adaptFromTask<DigitReader>(useMC, verbose)},

    Options{{"hmpid-digit-infile" /*"/qc-hmpid-digits"*/, VariantType::String, "hmpiddigits.root", {"Name of the input file with digits"}},
            {"input-dir", VariantType::String, "./", {"Input directory"}}}};
}

} // namespace hmpid

} // namespace o2
