// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   digits-to-cluster-workflow.cxx
/// \author Antonio Franco - INFN Bari
/// \version 1.0
/// \date 22 nov 2021
///

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/NameConf.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/CallbackService.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/ControlService.h"
#include "Framework/DataSpecUtils.h"
#include "Framework/DispatchPolicy.h"
#include "Framework/Task.h"
#include "Framework/Variant.h"
#include "Framework/WorkflowSpec.h"

/*
 ef : perform clusterization:
      either based on simulated data from a file, or on real data trhough a stream
      The executable reads upstream and writes upstream by default.
*/

void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}

// customize the completion policy
void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  using o2::framework::CompletionPolicy;
  using o2::framework::CompletionPolicyHelpers;
  policies.push_back(o2::framework::CompletionPolicyHelpers::defineByName(
    "clusters-hmpid-write", CompletionPolicy::CompletionOp::Consume));
}

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::string keyvaluehelp("Semicolon separated key=value strings ...");
  workflowOptions.push_back(
    o2::framework::ConfigParamSpec{"configKeyValues",
                                   o2::framework::VariantType::String,
                                   "",
                                   {keyvaluehelp}});
  workflowOptions.push_back(
    o2::framework::ConfigParamSpec{"read-from-file",
                                   o2::framework::VariantType::Bool,
                                   false,
                                   {"read upstream by default"}});
  workflowOptions.push_back(
    o2::framework::ConfigParamSpec{"write-to-file",
                                   o2::framework::VariantType::Bool,
                                   false,
                                   {"read upstream by default"}});

  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions);
}

#include "Framework/runDataProcessing.h"
#include "HMPIDWorkflow/DigitsToClustersSpec.h"
#include "HMPIDWorkflow/ClustersWriterSpec.h"
//#include "HMPIDWorkflow/HMPIDDigitizerSpec.h"

using namespace o2;
using namespace o2::framework;

WorkflowSpec defineDataProcessing(const ConfigContext& configcontext)
{
  WorkflowSpec specs;
  o2::conf::ConfigurableParam::updateFromString(
    configcontext.options().get<std::string>("configKeyValues"));
  auto mFromFile = configcontext.options().get<bool>(
    "read-from-file"); // read upstream by default
  auto mToFile = configcontext.options().get<bool>(
    "write-to-file"); // write upstream by default

  DataProcessorSpec consumer =
    o2::hmpid::getDigitsToClustersSpec("HMP/DIGITS", mFromFile, mToFile);
  specs.push_back(consumer);

  if (mToFile) { // Write to File
    DataProcessorSpec consumerClusterToRoot = o2::hmpid::getClustersToRootWriter();
    specs.push_back(consumerClusterToRoot);
  }

  return specs;
}
