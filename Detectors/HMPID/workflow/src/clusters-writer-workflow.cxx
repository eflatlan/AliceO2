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
/// \brief Workflow for clusterization for HMPID; read upstream/from file write upstream/to file.

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

void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  // ordered policies for the writers
  policies.push_back(o2::framework::CompletionPolicyHelpers::consumeWhenAllOrdered(".*HMPClustersWriter.*"));
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

  // ef : added
  workflowOptions.push_back(
    o2::framework::ConfigParamSpec{"disable-mc",
                                   o2::framework::VariantType::Bool,
                                   false,
                                   {"Do not propagate MC info"}});


  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions);
}

#include "Framework/runDataProcessing.h"
#include "HMPIDWorkflow/ClustersWriterSpec.h"

using namespace o2;
using namespace o2::framework;

WorkflowSpec defineDataProcessing(const ConfigContext& configcontext)
{
  WorkflowSpec specs;
  o2::conf::ConfigurableParam::updateFromString(
    configcontext.options().get<std::string>("configKeyValues"));

  // ef : added 
  bool useMC = !configcontext.options().get<bool>("disable-mc");
  specs.push_back(hmpid::getClusterWriterSpec(useMC));

  return specs;
}
