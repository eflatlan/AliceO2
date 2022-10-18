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

/// \file   DigitsReaderSpec.cxx
/// \author Annalisa Mastroserio - INFN Bari
/// \version 1.0
/// \date 22 Jun 2022
/// \brief Implementation of a data processor to read Digits tree and provide the array for further usage
///

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
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/Logger.h"
#include "Framework/DataRefUtils.h"
#include "Framework/InputRecordWalker.h"

#include "Headers/RAWDataHeader.h"
#include "DetectorsRaw/RDHUtils.h"
#include "DPLUtils/DPLRawParser.h"

//#include "DataFormatsHMP/Trigger.h"
//#include "DataFormatsHMP/Digits.h"
#include "HMPIDBase/Geo.h"

#include "HMPIDWorkflow/DigitsReaderSpec.h"

namespace o2
{
namespace hmpid
{

using namespace o2;
/// using namespace o2::header;
using namespace o2::framework;
using RDH = o2::header::RDHAny;

//
void DigitsReaderTask::init(framework::InitContext& ic)
{
  LOG(info) << "[HMPID Digits reader - init() ] ";
  mDigitsReceived = 0;
  // Build the file name
  const auto filename = o2::utils::Str::concat_string(
    o2::utils::Str::rectifyDirectory(
      ic.options().get<std::string>("input-dir")),
    ic.options().get<std::string>("qc-hmpid-digits"));
  initFileIn(filename);
}

// return;

void DigitsReaderTask::run(framework::ProcessingContext& pc)
{
  // check if more entries in tree

  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
    pc.services().get<ControlService>().endOfStream();
    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    mExTimer.stop();
    mExTimer.logMes("End DigitsReader !  digits = " +
                    std::to_string(mDigitsReceived));
  } else {
    auto entry = mTree->GetReadEntry() + 1;
    assert(entry < mTree->GetEntries());
    mTree->GetEntry(entry);
    pc.outputs().snapshot(Output{"HMP", "DIGITS", 0, Lifetime::Timeframe}, mDigitsFromFile);
    pc.outputs().snapshot(Output{"HMP", "INTRECORDS", 0, Lifetime::Timeframe}, mDigitsTriggerFromFile);
    mClustersReceived += mDigitsFromFile.size();
    LOG(info) << "[HMPID DigitsReader - run() ] digits  = " << mDigitsFromFile.size();
  }

  // mExTimer.elapseMes("# received Clusters = " + std::to_string(mClustersReceived));
  return;
}

/*svoid ClusterReaderTask::endOfStream(framework::EndOfStreamContext& ec)
{
 mExTimer.stop();
 mExTimer.logMes("End ClusterReader !  clusters = " +
                  std::to_string(mClustersReceived));
 return;
} */

void DigitsReaderTask::initFileIn(const std::string& filename)
{
  // Create the TFIle
  mTree.reset(nullptr);
  mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());
  mTree.reset((TTree*)mFile->Get("o2sim"));

  if (!mTree) {
    LOG(error)
      << "HMPID Digitstr();
    throw std::runtime_error(
      "HMPID DigitsReaderTask::init() : Did not find "
      "o2sim file in clusters tree");
  }

  mTree->SetBranchAddress("HMPIDdigits", &mDigitsFromFilePtr);
  mTree->SetBranchAddress("InteractionRecords", &mDigitsTriggersFromFilePtr);
  mTree->Print("toponly");
}

//_________________________________________________________________________________________________

o2::framework::DataProcessorSpec getDigitsReaderSpec()
{

  std::vector<o2::framework::OutputSpec> outputs;
  outputs.emplace_back("HMP", "DIGITS", 0, o2::framework::Lifetime::Timeframe);
  outputs.emplace_back("HMP", "INTRECORDS", 0, o2::framework::Lifetime::Timeframe);

  return DataProcessorSpec{
    "HMP-DigitReader",
    Inputs{},
    outputs,
    AlgorithmSpec{adaptFromTask<ClusterReaderTask>()},
    Options{{"qc-hmpid-digits", VariantType::String, "hmpiddigits.root", {"Name of the input file with digits"}},
            {"input-dir", VariantType::String, "./", {"Input directory"}}}};
}

} // namespace hmpid
} // end namespace o2
