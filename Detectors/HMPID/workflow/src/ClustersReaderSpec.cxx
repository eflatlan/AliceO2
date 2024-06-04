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


  //ef : for debugging, check how many different evendIDs for the same HMPID trigger
  if (mVerbose) {
    int tNum = 0;
    for (const auto trig : *mClusterTriggersFromFilePtr) {

      auto timeA = o2::InteractionRecord::bc2ns(trig.getBc(), trig.getOrbit());
      int cnt = 0;

      const int firstEntry = trig.getFirstEntry();
      const int lastEntry = trig.getLastEntry();



      if (mUseMC) {
        std::vector<int> eventLabels;

        for (int i = firstEntry; i <= lastEntry; i++) {

          if (i < mLabels.getIndexedSize() && i < mClustersFromFile.size()) {
            const auto& labels = mLabels.getLabels(i);
            int prevEventLabel = 0;

            if (labels.size() > 0) {
              prevEventLabel = labels[0].getEventID();
              eventLabels.push_back(prevEventLabel);
            }

            for (const auto& label : labels) {

              if (label.getEventID() != prevEventLabel) {
                eventLabels.push_back(label.getEventID());
              }

              prevEventLabel = label.getEventID();
            }
          }
        }



        LOGP(info, "Trigger number {} : entries {}", tNum, trig.getNumberOfObjects());
        LOGP(info, "\n Different labels from eventLabels {} :::", eventLabels.size());

        std::vector<int> sortedVec = eventLabels;
        std::sort(sortedVec.begin(), sortedVec.end());

        std::cout << "eventLabels values: ";
        for (size_t i = 0; i < sortedVec.size(); ++i) {
          if (i == sortedVec.size() - 1 || sortedVec[i] != sortedVec[i + 1]) {
            std::cout << sortedVec[i] << " , ";
          }
        }
      }
      tNum++;
    }
  }

  if (mUseMC) {
    pc.outputs().snapshot(Output{"HMP", "CLUSTERSMCTR", 0}, mLabels);
  }

  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
    pc.services().get<ControlService>().endOfStream();
    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    mExTimer.stop();
    mExTimer.logMes("End ClusterReader !  clusters = " +
                    std::to_string(mClustersReceived));
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

  // ef: get useMC, adpted from CPV
  if (mUseMC) {

    if (mTree->GetBranch(mClusterMCTruthBranchName.c_str())) {
      mTree->SetBranchAddress(mClusterMCTruthBranchName.c_str(), &mLabelsPtr);
    } else {
      LOG(warning) << "MC-truth is missing, message will be empty";
    }
  }

  mTree->Print("toponly");
  mTree->SetBranchAddress("InteractionRecords", &mClusterTriggersFromFilePtr);
}

//_________________________________________________________________________________________________

o2::framework::DataProcessorSpec getClusterReaderSpec(bool useMC, bool verbose)
{

  std::vector<o2::framework::OutputSpec> outputs;
  outputs.emplace_back("HMP", "CLUSTERS", 0, o2::framework::Lifetime::Timeframe);
  outputs.emplace_back("HMP", "INTRECORDS1", 0, o2::framework::Lifetime::Timeframe);

  // ef: added here
  if (useMC) {
    outputs.emplace_back("HMP", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "HMP-ClusterReader",
    Inputs{},
    outputs,
    AlgorithmSpec{adaptFromTask<ClusterReaderTask>(useMC, verbose)},
    Options{{"hmpid-cluster-infile" /*"qc-hmpid-clusters"*/, VariantType::String, "hmpidclusters.root", {"Name of the input file with clusters"}},
            {"input-dir", VariantType::String, "./", {"Input directory"}}}};
}

} // namespace hmpid
} // end namespace o2
