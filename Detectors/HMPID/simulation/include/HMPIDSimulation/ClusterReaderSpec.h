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

///
/// \file    DatDecoderSpec.h
/// \author  Andrea Ferrero
///
/// \brief Definition of a data processor to run the raw decoding
///

#ifndef DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_CLUSTERREADERSPEC_H_
#define DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_CLUSTERREADERSPEC_H_

// ROOT
#include <TFile.h>
#include <TTree.h>

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"

#include "HMPIDBase/Common.h"
#include "HMPIDReconstruction/Clusterer.h"
#include <DataFormatsHMP/Cluster.h>
#include <DataFormatsHMP/Trigger.h>
#include "Framework/WorkflowSpec.h"

namespace o2
{
namespace hmpid
{

class ClusterReaderTask : public framework::Task
{
 public:
  ClusterReaderTask() = default;
  ~ClusterReaderTask() override = default;

  void init(framework::InitContext& ic) final;

  void run(framework::ProcessingContext& pc) final;
  void endOfStream(framework::EndOfStreamContext& ec) override;

 private:
  void initFileIn(const std::string& filename) final;
  // std::string mSigmaCutPar;
  // float mSigmaCut[7] = {4, 4, 4, 4, 4, 4, 4};

  // o2::hmpid::Clusterer* mRec;
  // long mDigitsReceived;

  ExecutionTimer mExTimer;

  std::unique_ptr<TFile> mFile = nullptr; // root file with Clusters
  TTree* mTree = nullptr;                 // tree inside the file
  //std::vector<o2::hmpid::Trigger> intrecords; //, *pintrecords= &intrecords; // pointer to InteractionRecords branch
  //std::vector<o2::hmpid::Cluster> clusters; //, *pclusters = &clusters; // pointer to HMPIDCluster branch

  unsigned long mNumberOfEntries = 0; // number of entries from TTree
  unsigned long mCurrentEntry = 0;    // index of current entry

  //void strToFloatsSplit(std::string s, std::string delimiter, float* res, int maxElem = 7);
};

o2::framework::DataProcessorSpec getClusterReaderSpec();

} // end namespace hmpid
} // end namespace o2

#endif
