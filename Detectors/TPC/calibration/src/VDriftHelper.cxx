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

#include "TPCBase/CDBInterface.h"
#include "TPCCalibration/VDriftHelper.h"
#include "DataFormatsTPC/LtrCalibData.h"
#include "TPCBase/ParameterGas.h"
#include "TPCBase/ParameterDetector.h"
#include "TPCBase/ParameterElectronics.h"
#include "Framework/Logger.h"
#include "Framework/ProcessingContext.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/InputRecord.h"
#include "Framework/ConcreteDataMatcher.h"

using namespace o2::tpc;
using namespace o2::framework;

//________________________________________________________
VDriftHelper::VDriftHelper()
{
  const auto& gaspar = o2::tpc::ParameterGas::Instance();
  const auto& detpar = o2::tpc::ParameterDetector::Instance();
  const auto& elpar = o2::tpc::ParameterElectronics::Instance();
  mVD.corrFact = 1.0;
  mVD.refVDrift = gaspar.DriftV;
  mVD.refTimeOffset = detpar.DriftTimeOffset * elpar.ZbinWidth; // convert time bins to \mus
  // was it imposed from the command line?
  mVD.creationTime = 1;                                                                                                         // just to be above 0
  if (o2::conf::ConfigurableParam::getProvenance("TPCGasParam.DriftV") == o2::conf::ConfigurableParam::EParamProvenance::kRT) { // we stick to this value
    mVD.creationTime = std::numeric_limits<long>::max();
    mForceParamDrift = true;
    LOGP(info, "TPC VDrift was set from command line to {}, will neglect update from CCDB", mVD.refVDrift);
  }
  if (o2::conf::ConfigurableParam::getProvenance("TPCDetParam.DriftTimeOffset") == o2::conf::ConfigurableParam::EParamProvenance::kRT) { // we stick to this value
    mVD.creationTime = std::numeric_limits<long>::max();
    mForceParamOffset = true;
    LOGP(info, "TPC dridt time offset was set from command line to {} mus ({} TB), will neglect update from CCDB",
         mVD.refTimeOffset, detpar.DriftTimeOffset);
  }

  mUpdated = true;
  mSource = Source::Param;
}

//________________________________________________________
void VDriftHelper::accountLaserCalibration(const LtrCalibData* calib, long fallBackTimeStamp)
{
  if (!calib || mForceParamDrift) { // laser may set only DriftParam (the offset is 0)
    return;
  }
  // old entries of laser calib have no update time assigned
  long updateTS = calib->creationTime > 0 ? calib->creationTime : fallBackTimeStamp;
  LOG(info) << "accountLaserCalibration " << calib->refVDrift << " / " << calib->getDriftVCorrection() << " t " << updateTS << " vs " << mVDLaser.creationTime;
  // old entries of laser calib have no reference assigned
  float ref = calib->refVDrift > 0. ? calib->refVDrift : o2::tpc::ParameterGas::Instance().DriftV;
  float corr = calib->getDriftVCorrection();
  if (corr > 0.) { // laser correction is inverse multiplicative
    static bool firstCall = true;
    auto prevRef = mVDLaser.refVDrift;
    mVDLaser.refVDrift = ref;
    mVDLaser.corrFact = 1. / corr;
    mVDLaser.creationTime = calib->creationTime;
    mUpdated = true;
    mSource = Source::Laser;
    if (mMayRenormSrc & (0x1U << Source::Laser)) { // this was 1st setting?
      if (corr != 1.f) {                           // this may happen if old-style (non-normalized) standalone or non-normalized run-time laset calibration is used
        LOGP(warn, "VDriftHelper: renorming initinal TPC refVDrift={}/correction={} to {}/1.0, source: {}", mVDLaser.refVDrift, mVDLaser.corrFact, mVDLaser.getVDrift(), getSourceName(mSource));
        mVDLaser.normalize(); // renorm reference to have correction = 1.
      }
      mMayRenormSrc &= ~(0x1U << Source::Laser); // unset MayRenorm
    } else if (ref != prevRef) {                 // we want to keep the same reference over the run, this may happen if run-time laser calibration is supplied
      LOGP(warn, "VDriftHelper: renorming updated TPC refVDrift={}/correction={} previous refVDrift {}, source: {}", mVDLaser.refVDrift, mVDLaser.corrFact, prevRef, getSourceName(mSource));
      mVDLaser.normalizeOffset(prevRef);
    }
  }
}

//________________________________________________________
void VDriftHelper::accountDriftCorrectionITSTPCTgl(const VDriftCorrFact* calib)
{
  if (!calib || (mForceParamDrift && mForceParamOffset)) {
    return;
  }
  LOG(info) << "accountDriftCorrectionITSTPCTgl " << calib->corrFact << " t " << calib->creationTime << " vs " << mVDTPCITSTgl.creationTime;
  auto prevRefVDrift = mVDTPCITSTgl.refVDrift;
  auto prevRefTOffs = mVDTPCITSTgl.refTimeOffset;
  mVDTPCITSTgl = *calib;
  mUpdated = true;
  mSource = Source::ITSTPCTgl;
  if (mMayRenormSrc & (0x1U << Source::ITSTPCTgl)) {         // this was 1st setting?
    if (!mForceParamDrift && mVDTPCITSTgl.corrFact != 1.f) { // this may happen if calibration from prevous run is used
      LOGP(warn, "VDriftHelper: renorming initinal TPC refVDrift={}/correction={} to {}/1.0, source: {}", mVDTPCITSTgl.refVDrift, mVDTPCITSTgl.corrFact, mVDTPCITSTgl.getVDrift(), getSourceName(mSource));
      mVDTPCITSTgl.normalize(); // renorm reference to have correction = 1.
    }
    if (!mForceParamOffset && mVDTPCITSTgl.timeOffsetCorr != 0.) {
      LOGP(warn, "VDriftHelper: renorming initinal TPC refTimeOffset={}/correction={} to {}/0.0, source: {}", mVDTPCITSTgl.refTimeOffset, mVDTPCITSTgl.timeOffsetCorr, mVDTPCITSTgl.getTimeOffset(), getSourceName());
      mVDTPCITSTgl.normalizeOffset();
    }
    mMayRenormSrc &= ~(0x1U << Source::ITSTPCTgl); // unset MayRenorm
  } else {
    if (!mForceParamDrift && mVDTPCITSTgl.refVDrift != prevRefVDrift) { // we want to keep the same reference over the run, this should not happen!
      LOGP(warn, "VDriftHelper: renorming updated TPC refVDrift={}/correction={} previous refVDrift {}, source: {}", mVDTPCITSTgl.refVDrift, mVDTPCITSTgl.corrFact, prevRefVDrift, getSourceName());
      mVDTPCITSTgl.normalize(prevRefVDrift);
    }
    if (!mForceParamOffset && mVDTPCITSTgl.refTimeOffset != prevRefTOffs) { // we want to keep the same reference over the run, this should not happen!
      LOGP(warn, "VDriftHelper: renorming updated TPC refTimeOffset={}/correction={} previous refTimeOffset {}, source: {}", mVDTPCITSTgl.refTimeOffset, mVDTPCITSTgl.timeOffsetCorr, prevRefTOffs, getSourceName());
      mVDTPCITSTgl.normalizeOffset(prevRefTOffs);
    }
  }
}

//________________________________________________________
void VDriftHelper::extractCCDBInputs(ProcessingContext& pc, bool laser, bool itstpcTgl)
{
  if (mForceParamDrift && mForceParamOffset) { // fixed from the command line
    return;
  }
  if (laser && !mForceParamDrift) {
    pc.inputs().get<o2::tpc::LtrCalibData*>("laserCalib");
  }
  if (itstpcTgl) {
    pc.inputs().get<o2::tpc::VDriftCorrFact*>("vdriftTgl");
  }
  if (mUpdated) { // there was a change
    // prefer among laser and tgl VDrift the one with the latest update time
    auto saveVD = mVD;
    mVD = mVDTPCITSTgl.creationTime < mVDLaser.creationTime ? mVDLaser : mVDTPCITSTgl;
    if (mForceParamDrift) {
      mVD.refVDrift = saveVD.refVDrift;
      mVD.corrFact = saveVD.corrFact;
      mVD.corrFactErr = 0.f;
    }
    if (mForceParamOffset) {
      mVD.refTimeOffset = saveVD.refTimeOffset;
      mVD.timeOffsetCorr = 0.f;
    }
    mSource = mVDTPCITSTgl.creationTime < mVDLaser.creationTime ? Source::Laser : Source::ITSTPCTgl;
    LOGP(info, "Will prefer TPC Drift from {} with time {} to {} with time {}",
         SourceNames[int(mSource)], mVD.creationTime,
         mSource == Source::Laser ? SourceNames[int(Source::ITSTPCTgl)] : SourceNames[int(Source::Laser)],
         mSource == Source::Laser ? mVDTPCITSTgl.creationTime : mVDLaser.creationTime);
    if (mForceParamDrift || mForceParamOffset) {
      LOGP(info, "but {} is imposed from the command line", mForceParamDrift ? "VDrift" : "DriftTimeOffset");
    }
  }
}

//________________________________________________________
void VDriftHelper::requestCCDBInputs(std::vector<InputSpec>& inputs, bool laser, bool itstpcTgl)
{
  if (laser) {
    addInput(inputs, {"laserCalib", "TPC", "CalibLaserTracks", 0, Lifetime::Condition, ccdbParamSpec(CDBTypeMap.at(CDBType::CalLaserTracks))});
  }
  if (itstpcTgl) {
    // VDrift calibration may change during the run (in opposite to Laser calibration, at least at the moment), so ask per-TF query
    addInput(inputs, {"vdriftTgl", "TPC", "VDriftTgl", 0, Lifetime::Condition, ccdbParamSpec(CDBTypeMap.at(CDBType::CalVDriftTgl), {}, 1)});
  }
}

//________________________________________________________
void VDriftHelper::addInput(std::vector<InputSpec>& inputs, InputSpec&& isp)
{
  if (std::find(inputs.begin(), inputs.end(), isp) == inputs.end()) {
    inputs.emplace_back(isp);
  }
}

//________________________________________________________
bool VDriftHelper::accountCCDBInputs(const ConcreteDataMatcher& matcher, void* obj)
{
  if (matcher == ConcreteDataMatcher("TPC", "VDriftTgl", 0)) {
    accountDriftCorrectionITSTPCTgl(static_cast<VDriftCorrFact*>(obj));
    return true;
  }
  if (matcher == ConcreteDataMatcher("TPC", "CalibLaserTracks", 0)) {
    accountLaserCalibration(static_cast<LtrCalibData*>(obj));
    return true;
  }
  return false;
}
