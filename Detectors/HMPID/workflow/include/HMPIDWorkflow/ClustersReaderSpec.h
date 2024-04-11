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
/// \brief Definition of a data processor to read Clusters
///

#ifndef DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_CLUSTERSREADERSPEC_H_
#define DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_CLUSTERSREADERSPEC_H_

// ROOT
#include <TFile.h>
#include <TTree.h>

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "Framework/WorkflowSpec.h"

#include "HMPIDBase/Common.h"
#include "HMPIDBase/Geo.h"
#include "HMPIDReconstruction/Clusterer.h"
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Trigger.h"

namespace o2
{
namespace hmpid
{

class ClusterReaderTask : public framework::Task
{
 public:

  // ef added
  ClusterReaderTask(bool useMC, bool verbose)
  {
    mUseMC = useMC;
    mVerbose = verbose;
  };

  ~ClusterReaderTask() override = default;

  void init(framework::InitContext& ic) final;

  void run(framework::ProcessingContext& pc) final;
  // void endOfStream(framework::EndOfStreamContext& ec) override;

 private:
  // ef : added
  bool mUseMC = false;
  bool mVerbose = false;

  std::string mClusterMCTruthBranchName = "HMPIDClusterLabels";

  bool mReadFile = false;
  void initFileIn(const std::string& filename);

  long mClustersReceived;
  ExecutionTimer mExTimer;

  std::unique_ptr<TFile> mFile; // root file with Clusters
  std::unique_ptr<TTree> mTree; // tree inside the file

  // ef: add mLabels for clusteres
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mLabels, *mLabelsPtr = &mLabels;

  std::vector<o2::hmpid::Trigger> mClusterTriggersFromFile, *mClusterTriggersFromFilePtr = &mClusterTriggersFromFile;
  std::vector<o2::hmpid::Cluster> mClustersFromFile, *mClustersFromFilePtr = &mClustersFromFile;

  unsigned long mNumberOfEntries = 0; // number of entries from TTree
  unsigned long mCurrentEntry = 0;    // index of current entry

  // void strToFloatsSplit(std::string s, std::string delimiter, float* res, int maxElem = 7);
};

o2::framework::DataProcessorSpec getClusterReaderSpec(bool useMC, bool verbose = false);

} // end namespace hmpid
} // end namespace o2

#endif
