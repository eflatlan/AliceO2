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

/// @file   HMPMatcherSpec.cxx

#include <vector>
#include <string>
#include "TStopwatch.h"
#include "Framework/ConfigParamRegistry.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "Framework/Task.h"
#include "Framework/DataProcessorSpec.h"

// from Tracks
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/GlobalTrackAccessor.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsITS/TrackITS.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "ReconstructionDataFormats/TrackTPCTOF.h"

// from HMPID
#include "DataFormatsHMP/Cluster.h"
#include "GlobalTracking/MatchHMP.h"
#include "GlobalTrackingWorkflow/HMPMatcherSpec.h"

using namespace o2::framework;
// using MCLabelsTr = gsl::span<const o2::MCCompLabel>;
// using GID = o2::dataformats::GlobalTrackID;
// using DetID = o2::detectors::DetID;

using evIdx = o2::dataformats::EvIndex<int, int>;
using MatchOutputType = std::vector<o2::dataformats::MatchInfoHMP>;
using GID = o2::dataformats::GlobalTrackID;

namespace o2
{
namespace globaltracking
{

class HMPMatcherSpec : public Task
{
 public:
  HMPMatcherSpec(std::shared_ptr<DataRequest> dr, bool useMC, bool useFIT, bool tpcRefit, bool strict) : mDataRequest(dr), mUseMC(useMC), mUseFIT(useFIT), mDoTPCRefit(tpcRefit), mStrict(strict) {}
  ~HMPMatcherSpec() override = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(framework::EndOfStreamContext& ec) final;

 private:
  std::shared_ptr<DataRequest> mDataRequest;
  bool mUseMC = true;
  bool mUseFIT = false;
  bool mDoTPCRefit = false;
  bool mStrict = false;
  MatchHMP mMatcher; ///< Cluster finder
  TStopwatch mTimer;
};

void HMPMatcherSpec::init(InitContext& ic)
{
  mTimer.Stop();
  mTimer.Reset();
  //-------- init geometry and field --------//
  o2::base::GeometryManager::loadGeometry();
  o2::base::Propagator::initFieldFromGRP();
  std::unique_ptr<o2::parameters::GRPObject> grp{o2::parameters::GRPObject::loadFrom()};

  // this is a hack to provide Mat.LUT from the local file, in general will be provided by the framework from CCDB
  std::string matLUTPath = ic.options().get<std::string>("material-lut-path");
  std::string matLUTFile = o2::base::NameConf::getMatLUTFileName(matLUTPath);
  if (o2::utils::Str::pathExists(matLUTFile)) {
    auto* lut = o2::base::MatLayerCylSet::loadFromFile(matLUTFile);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    LOG(debug) << "Loaded material LUT from " << matLUTFile;
  } else {
    LOG(debug) << "Material LUT " << matLUTFile << " file is absent, only TGeo can be used";
  }
  if (mStrict) {
    // mMatcher.setHighPurity(); ef: function is commented out in MatchHMP.h
  }
}

void HMPMatcherSpec::run(ProcessingContext& pc)
{
  mTimer.Start(false);

  RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());

  auto creationTime = DataRefUtils::getHeader<DataProcessingHeader*>(pc.inputs().getFirstValid(true))->creation;

  LOG(debug) << "isTrackSourceLoaded: TPC -> " << recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::TPC);
  LOG(debug) << "isTrackSourceLoaded: ITSTPC -> " << recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::ITSTPC);
  LOG(debug) << "isTrackSourceLoaded: TPCTRD -> " << recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::TPCTRD);
  LOG(debug) << "isTrackSourceLoaded: ITSTPCTRD -> " << recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::ITSTPCTRD);

  bool isTPCused = recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::TPC);
  bool isITSTPCused = recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::ITSTPC);
  bool isTPCTRDused = recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::TPCTRD);
  bool isITSTPCTRDused = recoData.isTrackSourceLoaded(o2::dataformats::GlobalTrackID::Source::ITSTPCTRD);
  uint32_t ss = o2::globaltracking::getSubSpec(mStrict ? o2::globaltracking::MatchingType::Strict : o2::globaltracking::MatchingType::Standard);

  // mMatcher.setFIT(mUseFIT);

  // mMatcher.setTS(creationTime);

  mMatcher.run(recoData);

  if (isTPCused) {
    pc.outputs().snapshot(Output{o2::header::gDataOriginTOF, "MTC_TPC", ss, Lifetime::Timeframe}, mMatcher.getMatchedTrackVector(o2::dataformats::MatchInfoTOFReco::TrackType::TPC));
    if (mUseMC) {
      pc.outputs().snapshot(Output{o2::header::gDataOriginTOF, "MCMTC_TPC", ss, Lifetime::Timeframe}, mMatcher.getMatchedHMPLabelsVector(o2::dataformats::MatchInfoTOFReco::TrackType::TPC));
      // ef : what should the namescope of TrackType be?
    } // ef changed to getMatchedHMPLabelsVector

    auto nmatch = mMatcher.getMatchedTrackVector(o2::dataformats::MatchInfoTOFReco::TrackType::TPC).size();
    if (mDoTPCRefit) {
      LOG(debug) << "Refitting " << nmatch << " matched TPC tracks with TOF time info";
    } else {
      LOG(debug) << "Shifting Z for " << nmatch << " matched TPC tracks according to TOF time info";
    }
    auto& tracksTPCTOF = pc.outputs().make<std::vector<o2::dataformats::TrackTPCTOF>>(OutputRef{"tpctofTracks", ss}, nmatch);
    // mMatcher.makeConstrainedTPCTracks(tracksTPCTOF);
  }

  if (isITSTPCused) {
    pc.outputs().snapshot(Output{o2::header::gDataOriginHMP, "MTC_ITSTPC", 0, Lifetime::Timeframe}, mMatcher.getMatchedTrackVector(o2::dataformats::MatchInfoTOFReco::TrackType::ITSTPC));
    if (mUseMC) {
      pc.outputs().snapshot(Output{o2::header::gDataOriginHMP, "MCMTC_ITSTPC", 0, Lifetime::Timeframe}, mMatcher.getMatchedHMPLabelsVector(o2::dataformats::MatchInfoTOFReco::TrackType::ITSTPC));
    }
  }

  if (isTPCTRDused) {
    pc.outputs().snapshot(Output{o2::header::gDataOriginTOF, "MTC_TPCTRD", ss, Lifetime::Timeframe}, mMatcher.getMatchedTrackVector(o2::dataformats::MatchInfoTOFReco::TrackType::TPCTRD));
    if (mUseMC) {
      pc.outputs().snapshot(Output{o2::header::gDataOriginTOF, "MCMTC_TPCTRD", ss, Lifetime::Timeframe}, mMatcher.getMatchedHMPLabelsVector(o2::dataformats::MatchInfoTOFReco::TrackType::TPCTRD));
    }
  }

  if (isITSTPCTRDused) {
    pc.outputs().snapshot(Output{o2::header::gDataOriginTOF, "MTC_ITSTPCTRD", 0, Lifetime::Timeframe}, mMatcher.getMatchedTrackVector(o2::dataformats::MatchInfoTOFReco::TrackType::ITSTPCTRD));
    if (mUseMC) {
      pc.outputs().snapshot(Output{o2::header::gDataOriginTOF, "MCMTC_ITSTPCTRD", 0, Lifetime::Timeframe}, mMatcher.getMatchedHMPLabelsVector(o2::dataformats::MatchInfoTOFReco::TrackType::ITSTPCTRD));
    }
  }

  // TODO: TRD-matched tracks
  //  pc.outputs().snapshot(Output{o2::header::gDataOriginTOF, "CALIBDATA", 0, Lifetime::Timeframe}, mMatcher.getCalibVector());

  mTimer.Stop();
}

void HMPMatcherSpec::endOfStream(EndOfStreamContext& ec)
{
  LOGF(debug, "HMP matching total timing: Cpu: %.3e Real: %.3e s in %d slots",
       mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
}

DataProcessorSpec getHMPMatcherSpec(GID::mask_t src, bool useMC, bool useFIT, bool tpcRefit, bool strict)
{
  uint32_t ss = o2::globaltracking::getSubSpec(strict ? o2::globaltracking::MatchingType::Strict : o2::globaltracking::MatchingType::Standard);
  auto dataRequest = std::make_shared<DataRequest>();
  if (strict) {
    dataRequest->setMatchingInputStrict();
  }
  dataRequest->requestTracks(src, useMC);
  dataRequest->requestClusters(GID::getSourceMask(GID::TOF), useMC);
  if (useFIT) {
    dataRequest->requestClusters(GID::getSourceMask(GID::FT0), false);
  }

  std::vector<OutputSpec> outputs;
  if (GID::includesSource(GID::TPC, src)) {
    outputs.emplace_back(o2::header::gDataOriginHMP, "MTC_TPC", ss, Lifetime::Timeframe);
    // outputs.emplace_back(OutputLabel{"tpctofTracks"}, o2::header::gDataOriginTOF, "TOFTRACKS_TPC", ss, Lifetime::Timeframe);
    if (useMC) {
      outputs.emplace_back(o2::header::gDataOriginHMP, "MCMTC_TPC", ss, Lifetime::Timeframe);
    }
  }
  if (GID::includesSource(GID::ITSTPC, src)) {
    outputs.emplace_back(o2::header::gDataOriginHMP, "MTC_ITSTPC", 0, Lifetime::Timeframe);
    if (useMC) {
      outputs.emplace_back(o2::header::gDataOriginHMP, "MCMTC_ITSTPC", 0, Lifetime::Timeframe);
    }
  }
  if (GID::includesSource(GID::ITSTPCTRD, src)) {
    outputs.emplace_back(o2::header::gDataOriginHMP, "MTC_ITSTPCTRD", 0, Lifetime::Timeframe);
    if (useMC) {
      outputs.emplace_back(o2::header::gDataOriginHMP, "MCMTC_ITSTPCTRD", 0, Lifetime::Timeframe);
    }
  }
  if (GID::includesSource(GID::TPCTRD, src)) {
    outputs.emplace_back(o2::header::gDataOriginHMP, "MTC_TPCTRD", ss, Lifetime::Timeframe);
    if (useMC) {
      outputs.emplace_back(o2::header::gDataOriginHMP, "MCMTC_TPCTRD", ss, Lifetime::Timeframe);
    }
  }
  // outputs.emplace_back(o2::header::gDataOriginTOF, "CALIBDATA", 0, Lifetime::Timeframe);

  return DataProcessorSpec{
    "hmp-matcher",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<HMPMatcherSpec>(dataRequest, useMC, useFIT, tpcRefit, strict)},
    Options{
      {"material-lut-path", VariantType::String, "", {"Path of the material LUT file"}}}};
}

} // namespace globaltracking
} // namespace o2