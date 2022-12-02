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

/// \file digits-to-clusters-workflow.h
/// \brief Workflow for clusterization for HMPID; read upstream/from file write upstream/to file

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


using o2::framework::ConfigParamSpec;
using o2::framework::VariantType;
// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  std::string keyvaluehelp("Semicolon separated key=value strings ...");
  workflowOptions.push_back(
    ConfigParamSpec{"configKeyValues",
                    VariantType::String,
                    "",
                    {keyvaluehelp}});

  workflowOptions.push_back(
   ConfigParamSpec{"read-from-file",
                   o2::framework::VariantType::Bool,
                   false,
                   {"read upstream by default"}});


  workflowOptions.push_back(
    ConfigParamSpec{"write-to-file",
                    VariantType::Bool,
                    false,
                    {"read upstream by default"}});

  workflowOptions.push_back(
    ConfigParamSpec{"out-file",
                    o2::framework::VariantType::String,
                    "hmpidclusters.root",
                    {"Output file"}});

  workflowOptions.push_back(
    ConfigParamSpec{"out-dir",
                    VariantType::String,
                    "",
                    {"Output directory"}});


  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions);
}

#include "Framework/runDataProcessing.h"
#include "HMPIDWorkflow/DigitsToClustersSpec.h"
#include "HMPIDWorkflow/ClustersWriterSpec.h"

//#include "HMPIDWorkflow/HMPIDDigitizerSpec.h"ss

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

  auto outFile = configcontext.options().get<std::string>(
    "out-file"); // specify output file
  auto outDir = configcontext.options().get<std::string>(
    "out-dir");  // specify output directory

  DataProcessorSpec consumer =
    o2::hmpid::getDigitsToClustersSpec("HMP/DIGITS", mFromFile, mToFile);

  specs.push_back(consumer);

  if (mToFile) { // Write to File
    //LOGP(warn, "Writing to file _ cstr  {}", t.c_str());
    //LOGP(warn, "Writing to file _   {}", t);
    //LOGP(warn, "Writing to dir  {}", outDir);
    //LOGP(warn, "Writing to file  {}", outFile);
    DataProcessorSpec consumerClusterToRoot = o2::hmpid::getClustersToRootWriter(outDir /*= ""*/, outFile /* = "hmpidclusters.root"*/);
    specs.push_back(consumerClusterToRoot);
  }

  return specs;
}
