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

#include "HMPIDWorkflow/ClustersWriterSpec.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "Framework/InputSpec.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "DataFormatsHMP/Cluster.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

/*ef:comment
#include <bits/stdc++.h>
#include<sys/stat.h>
#include<sys/types.h>
*/ 

namespace o2
{
namespace hmpid
{

template <typename T>
using BranchDefinition = framework::MakeRootTreeWriterSpec::BranchDefinition<T>;

o2::framework::DataProcessorSpec getClustersToRootWriter(std::string outDir, std::string outFile)
{
  using InputSpec = framework::InputSpec;
  using MakeRootTreeWriterSpec = framework::MakeRootTreeWriterSpec;

  /*
  LOGP(warn, "output loc and file = {}", outDir);
  LOGP(warn, "output loc and file = {}", outFile);

  LOGP(warn, "output loc and file c_str = {}", outDir.c_str());
  LOGP(warn, "output loc and file c_str = {}", outFile.c_str()); */ 

  //ef:comment 
  auto check = mkdir(outDir.c_str(), 0777);
  //LOGP(warn, "check dir = {check}");

  auto output = "./" + outDir + "/" + outFile;
  LOGP(warn, "output loc and file = {}", output);
  LOGP(warn, "output loc and file c_str = {}", output.c_str());

  return MakeRootTreeWriterSpec("HMPClustersWriter",
                                outFile.c_str(),
                                "o2sim",
                                1,
                BranchDefinition<std::vector<o2::hmpid::Cluster>>{InputSpec{"hmpclusterinput", "HMP", "CLUSTERS"}, "HMPIDClusters"},
                                BranchDefinition<std::vector<o2::hmpid::Trigger>>{InputSpec{"hmpinteractionrecords", "HMP", "INTRECORDS1"}, "InteractionRecords"})();
}

// ef: corrected misspelt branch-name HMPIDclusters --> HMPIDClusters

} // end namespace hmpid
} // end namespace o2
