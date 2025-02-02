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

#ifndef O2_TRD_KRCLUSTERERSPEC_H
#define O2_TRD_KRCLUSTERERSPEC_H

/// \file   KrClustererSpec.h
/// \brief Steers the TRD Krypton cluster finder
/// \author Ole Schmidt

// input TRD digits, TRD trigger records
// output Kr clusters

#include "TRDWorkflow/KrClustererSpec.h"
#include "TRDCalibration/KrClusterFinder.h"
#include "Framework/Task.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "TStopwatch.h"
#include <fairlogger/Logger.h>

using namespace o2::framework;

namespace o2
{
namespace trd
{
class TRDKrClustererDevice : public Task
{
 public:
  TRDKrClustererDevice() = default;
  ~TRDKrClustererDevice() override = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(framework::EndOfStreamContext& ec) final;

 private:
  o2::trd::KrClusterFinder mKrClFinder;
};

void TRDKrClustererDevice::init(InitContext& ic)
{
  mKrClFinder.init();
}

void TRDKrClustererDevice::run(ProcessingContext& pc)
{
  TStopwatch timer;

  const auto digits = pc.inputs().get<gsl::span<Digit>>("digits");
  const auto triggerRecords = pc.inputs().get<gsl::span<TriggerRecord>>("triggerRecords");

  mKrClFinder.reset();
  mKrClFinder.setInput(digits, triggerRecords);
  timer.Start();
  mKrClFinder.findClusters();
  timer.Stop();

  LOGP(info, "Found {} Kr clusters in {} input trigger records. Timing: CPU: {}, Real: {}",
       mKrClFinder.getKrClusters().size(), triggerRecords.size(), timer.CpuTime(), timer.RealTime());

  pc.outputs().snapshot(Output{o2::header::gDataOriginTRD, "KRCLUSTER", 0}, mKrClFinder.getKrClusters());
  pc.outputs().snapshot(Output{o2::header::gDataOriginTRD, "TRGKRCLS", 0}, mKrClFinder.getKrTrigRecs());
}

void TRDKrClustererDevice::endOfStream(EndOfStreamContext& ec)
{
  LOG(info) << "Done with the cluster finding (EoS received)";
}

framework::DataProcessorSpec getKrClustererSpec()
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("digits", ConcreteDataTypeMatcher{o2::header::gDataOriginTRD, "DIGITS"}, Lifetime::Timeframe);
  inputs.emplace_back("triggerRecords", ConcreteDataTypeMatcher{o2::header::gDataOriginTRD, "TRKTRGRD"}, Lifetime::Timeframe);
  std::vector<OutputSpec> outputs;
  outputs.emplace_back(o2::header::gDataOriginTRD, "KRCLUSTER", 0, Lifetime::Timeframe);
  outputs.emplace_back(o2::header::gDataOriginTRD, "TRGKRCLS", 0, Lifetime::Timeframe);

  return DataProcessorSpec{
    "trd-kr-clusterer",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<o2::trd::TRDKrClustererDevice>()},
    Options{}};
}

} // namespace trd

} // namespace o2

#endif // O2_TRD_KRCLUSTERERSPEC_H
