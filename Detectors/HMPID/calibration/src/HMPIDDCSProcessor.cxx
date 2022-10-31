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

#include "HMPIDCalibration/HMPIDDCSProcessor.h"
#include <TCanvas.h>

namespace o2::hmpid
{

// initialize map of DPIDs, function is called once in the
// HMPIDDCSDataProcessorSpec::init() function
void HMPIDDCSProcessor::init(const std::vector<DPID>& pids)
{
  for (const auto& it : pids) {
    mPids[it] = false;
  }
  arNmean.resize(43);
  arQthre.resize(42);
}

// process span of Datapoints, function is called repeatedly in the
// HMPIDDCSDataProcessorSpec::run() function
void HMPIDDCSProcessor::process(const gsl::span<const DPCOM> dps)
{

  // if there is no entries in span
  if (dps.size() == 0) {
    LOG(debug) << "Size = 0: ";
    return;
  }

  if (mVerbose) {
    LOG(debug) << "\n\n\nProcessing new DCS DP map\n-----------------";
  }

  if (!mFirstTimeSet) {
    mFirstTime = mStartValidity;
    mFirstTimeSet = true;
  }

  // itterate over span of datapoints
  for (const auto& dp : dps) {
    const auto& el = mPids.find(dp.id);

    // check if datapoint is within the defined map of DPs for HMPID
    if (el == mPids.end()) {
      LOG(info) << "DP " << dp.id << "Not found, will not be processed";
      continue;
    }
    // The aliasstring has been processed:
    mPids[dp.id] = true;

    // extract alias-string of DataPoints
    const std::string_view alias(dp.id.get_alias());
    const auto detectorId = alias.substr(0, 4);      // HMP_
    const auto transparencyId = alias.substr(0, 22); // HMP_TRANPLANT_MEASURE_

    // check if given dp is from HMPID
    if (transparencyId == TRANS_ID) {
      processTRANS(dp);
    } else if (detectorId == HMPID_ID) {
      processHMPID(dp);
    } else { // ef: changed to warn and "unkown DP" ==> "missing DP"
      LOG(warn) << "Missing data point: " << alias;
    }
  } // end for
}

// if the string of the dp contains the HMPID-specifier "HMP_",
// but not the Transparency-specifier
void HMPIDDCSProcessor::processHMPID(const DPCOM& dp)
{

  const std::string alias(dp.id.get_alias());
  const auto hmpidString = alias.substr(alias.length() - 8);

  if (hmpidString == TEMP_IN_ID) {
    if (mVerbose) {
      LOG(info) << "Temperature_in DP: " << alias;
    }
    fillTempIn(dp);
  } else if (hmpidString == TEMP_OUT_ID) {
    if (mVerbose) {
      LOG(info) << "Temperature_out DP: " << alias;
    }
    fillTempOut(dp);
  } else if (hmpidString == HV_ID) {
    if (mVerbose) {
      LOG(info) << "HV DP: " << alias;
    }
    fillHV(dp);
  } else if (hmpidString == ENV_PRESS_ID) {
    if (mVerbose) {
      LOG(info) << "Environment Pressure DP: " << alias;
    }
    fillEnvPressure(dp);
  } else if (hmpidString == CH_PRESS_ID) {
    if (mVerbose) {
      LOG(info) << "Chamber Pressure DP: " << alias;
    }
    fillChPressure(dp);
  } else { // ef: changed to warn and unkown DP ==> missing data DP
    LOG(warn) << "Missing data point: " << alias;
  }
}

// if the string of the dp contains the Transparency-specifier
// "HMP_TRANPLANT_MEASURE_"
void HMPIDDCSProcessor::processTRANS(const DPCOM& dp)
{
  const auto& dpid = dp.id;
  const std::string alias(dpid.get_alias());
  const auto transparencyString = alias.substr(alias.length() - 9);

  // Get index [0..29] of Transparency-measurement
  const int dig1 = aliasStringToInt(dpid, 22);
  const int dig2 = aliasStringToInt(dpid, 23);
  const int num = dig1 * 10 + dig2;

  if (dig1 == -1 || dig2 == -1) {
    LOG(warn) << "digits in string invalid" << dig1 << " " << dig2;
    return;
  }

  if (num < 0 || num > 29) {
    LOG(warn) << "num out of range " << num;
    return;
  }

  if (alias.substr(alias.length() - 10) == WAVE_LEN_ID) {
    waveLenVec[num].emplace_back(dp);
    if (mVerbose) {
      LOG(info) << "WAVE_LEN_ID DP: " << alias;
    }
  } else if (transparencyString == FREON_CELL_ID) {
    freonCellVec[num].emplace_back(dp);
    if (mVerbose) {
      LOG(info) << "FREON_CELL_ID DP: " << alias;
    }
  } else if (transparencyString == ARGON_CELL_ID) {
    if (mVerbose) {
      LOG(info) << "ARGON_CELL_ID DP: " << alias;
    }
    argonCellVec[num].emplace_back(dp);
  } else if (transparencyString == REF_ID) {
    if (alias.substr(alias.length() - 14) == ARGON_REF_ID) {
      if (mVerbose) {
        LOG(info) << "ARGON_REF_ID DP: " << alias;
      }
      argonRefVec[num].emplace_back(dp);
    } else if (alias.substr(alias.length() - 14) == FREON_REF_ID) {
      freonRefVec[num].emplace_back(dp);
      if (mVerbose) {
        LOG(info) << "FREON_REF_ID DP: " << alias;
      }
    } else {  // ef: remove mVerbose, change to warn, and "DP not found" ==> "Missing DP"
      LOG(warn) << "Missing Data point: " << alias;
    }
  } else {  // ef:change to warn and "DP not found" ==> "Missing DP"
    LOG(warn) << "Missing Data point: " << alias;
  }
}

// ======Fill
// DPCOM-entries===========================================================

// fill entries in environment pressure DPCOM-vector
void HMPIDDCSProcessor::fillEnvPressure(const DPCOM& dpcom)
{
  auto& dpid = dpcom.id;
  const auto& type = dpid.get_type();

  // check if datatype is as expected
  if (type == DeliveryType::DPVAL_DOUBLE) {
    dpVecEnv.emplace_back(dpcom);
  } else {
    LOG(warn) << "Invalid Datatype for Env Pressure";
    LOG(warn) << "Env Pressure string: " << dpid;
  }
}

// fill entries in chamber-pressure DPCOM-vector
void HMPIDDCSProcessor::fillChPressure(const DPCOM& dpcom)
{

  auto& dpid = dpcom.id;
  const auto& type = dpid.get_type();

  if (type == DeliveryType::DPVAL_DOUBLE) {

    // find chamber number:
    auto chNum = aliasStringToInt(dpid, indexChPr);
    if (chNum < 7 && chNum >= 0) {
      dpVecCh[chNum].emplace_back(dpcom);
    } else {
      LOG(warn) << "Chamber Number out of range for Pressure : " << chNum;
      const std::string inputString(dpid.get_alias());
      char stringPos = inputString[indexChPr];
      LOG(warn) << "Chamber Pressure string: " << inputString;
      LOG(warn) << "Chamber Pressure extracted num: " << stringPos;
    }
  } else {
    LOG(warn) << "Not correct datatype for Pressure : ";
  }
}

// HV in each chamber_section = 7*3 --> will result in Q_thre
void HMPIDDCSProcessor::fillHV(const DPCOM& dpcom)
{
  auto& dpid = dpcom.id;
  const auto& type = dpid.get_type();

  if (type == DeliveryType::DPVAL_DOUBLE) {
    const auto chNum = aliasStringToInt(dpid, indexChHv);
    const auto secNum = aliasStringToInt(dpid, indexSecHv);

    // ef: can remove this now, since the problem with the strings is fixed:
    /* 
    if (mVerbose) {
      LOGP(info, "HV ch:{} sec:{} val: {}", chNum, secNum, o2::dcs::getValue<double>(dpcom));
    } */
    if (chNum < 7 && chNum >= 0) {
      if (secNum < 6 && secNum >= 0) {
        dpVecHV[6 * chNum + secNum].emplace_back(dpcom);
      } else {
        LOG(warn) << "Sector Number out of range for HV : " << secNum;
      }
    } else {
      LOG(warn) << "Chamber Number out of range for HV : " << chNum;
    }
  } else {
    LOG(warn) << "Not correct datatype for HV DP";
  }
}

// Temp in (T1)  in each chamber_radiator = 7*3
void HMPIDDCSProcessor::fillTempIn(const DPCOM& dpcom)
{
  auto& dpid = dpcom.id;
  const auto& type = dpid.get_type();

  if (type == DeliveryType::DPVAL_DOUBLE) {
    auto chNum = aliasStringToInt(dpid, indexChTemp);
    auto radNum = aliasStringToInt(dpid, indexRadTemp);

    // verify chamber- and raiator-numbers
    if (chNum < 7 && chNum >= 0) {
      if (radNum < 3 && radNum >= 0) {
        dpVecTempIn[3 * chNum + radNum].emplace_back(dpcom);
      } else {
        LOG(warn) << "Radiator Number out of range for TempIn :" << radNum;
      }
    } else {
      LOG(warn) << "Chamber Number out of range for TempIn DP :" << chNum;
    }
  } else {
    LOG(warn) << "Not correct datatype for TempIn DP ";
  }
}

// Temp out (T2), in each chamber_radiator = 7*3
void HMPIDDCSProcessor::fillTempOut(const DPCOM& dpcom)
{
  auto& dpid = dpcom.id;
  const auto& type = dpid.get_type();

  if (type == DeliveryType::DPVAL_DOUBLE) {
    auto chNum = aliasStringToInt(dpid, indexChTemp);
    auto radNum = aliasStringToInt(dpid, indexRadTemp);

    // verify chamber- and raiator-numbers
    if (chNum < 7 && chNum >= 0) {
      if (radNum < 3 && radNum >= 0) {
        dpVecTempOut[3 * chNum + radNum].emplace_back(dpcom);

      } else {
        LOG(warn) << "Radiator Number out of range for TempOut DP : " << radNum;
      }
    } else {
      LOG(warn) << "Chamber Number out of range for TempOut DP : " << chNum;
    }
  } else {
    LOG(warn) << "Not correct datatype for TempOut DP ";
  }
}

//==== Calculate mean photon energy=============================================

double HMPIDDCSProcessor::procTrans()
{
  for (int i = 0; i < 30; i++) {

    photEn = calculatePhotonEnergy(i);

    if (photEn < o2::hmpid::Param::ePhotMin() ||
        photEn > o2::hmpid::Param::ePhotMax()) {
      LOG(warn) << "photon energy out of range" << photEn;
      continue; // if photon energy is out of range
    }

    // ===== evaluate phototube current for argon reference ============================
    refArgon = dpVector2Double(argonRefVec[i], "ARGONREF", i);
    if (refArgon == eMeanDefault) { 
      // ef: removed warn here, since it is done in the dpVector2double
      /*LOG(warn) << "refArgon == defaultEMean()"; */
      // ef: simply return eMeanDefault; replacing defaultEMean() function-call
      return eMeanDefault;
    }

    //===== evaluate phototube current for argon cell===================================
    cellArgon = dpVector2Double(argonCellVec[i], "ARGONCELL", i);
    if (cellArgon == eMeanDefault) {
      return eMeanDefault;
    }

    //==== evaluate phototube current for freon reference ==============================
    refFreon = dpVector2Double(freonRefVec[i], "C6F14REFERENCE", i);
    if (refFreon == eMeanDefault) {
      return eMeanDefault;
    }

    // ==== evaluate phototube current for freon cell ==================================
    cellFreon = dpVector2Double(freonCellVec[i], "C6F14CELL", i);
    if (cellFreon == eMeanDefault) {
      return eMeanDefault;
    }

    // evaluate correction factor to calculate trasparency
    bool isEvalCorrOk = evalCorrFactor(refArgon, cellArgon, refFreon, cellFreon, photEn, i);

    // ef: Returns false if dRefFreon * dRefArgon < 0
    if (!isEvalCorrOk) {
      return eMeanDefault;
    }


    // Evaluate timestamps :

    auto s1 = getMinTime(waveLenVec[i]);
    auto s2 = getMinTime(argonRefVec[i]);
    auto s3 = getMinTime(argonCellVec[i]);
    auto s4 = getMinTime(freonRefVec[i]);
    auto s5 = getMinTime(freonCellVec[i]);

    auto e1 = getMaxTime(waveLenVec[i]);
    auto e2 = getMaxTime(argonRefVec[i]);
    auto e3 = getMaxTime(argonCellVec[i]);
    auto e4 = getMaxTime(freonRefVec[i]);
    auto e5 = getMaxTime(freonCellVec[i]);

    auto minTime = std::max({s1, s2, s3, s4, s5});
    auto maxTime = std::min({e1, e2, e3, e4, e5});

    if (minTime < mTimeEMean.first) {
      mTimeEMean.first = minTime;
    }
    if (maxTime > mTimeEMean.last) {
      mTimeEMean.last = maxTime;
    }

  } // end for


  // evaluate total energy --> mean photon energy
  if (sProb > 0) {
    eMean = sEnergProb / sProb;
  } else {
    // ef: removed mVerbose and changed to warn
    LOGP(warn, " sProb < 0 --> Default E mean used! ");
    return eMeanDefault;
  }

  // ef: changed to warn
  if (eMean < o2::hmpid::Param::ePhotMin() ||
      eMean > o2::hmpid::Param::ePhotMax()) { 
    LOGP(warn, "eMean out of range  ({}) --> Default E mean used! ", eMean);
    return eMeanDefault;
  }

  return eMean;

} // end ProcTrans

// ==== procTrans help-functions =============================

//ef: do not use function-call; just return the value instead
/*
double HMPIDDCSProcessor::defaultEMean()
{
  //ef: changed to warn, 
  //    not removing mVerbose, as it could be called very often
  if (mVerbose) {
    LOG(warn) << Form(" Mean energy photon calculated ---> %f eV ", eMeanDefault);
  }
  return eMeanDefault;
} */

//==== evaluate wavelenght
//=======================================================
double HMPIDDCSProcessor::calculatePhotonEnergy(int i)
{
  // if there is no entries
  if (waveLenVec[i].size() == 0) {
    LOG(warn) << Form("No Data Point values for %i.waveLenght --> Default E mean used!", i);
    return eMeanDefault; // will break this entry in foor loop
  }

  DPCOM dp = (waveLenVec[i])[0];

  // check if datatype is as expected
  if (dp.id.get_type() == DeliveryType::DPVAL_DOUBLE) {
    lambda = o2::dcs::getValue<double>(dp);
  } else {
    LOG(warn) << Form("Not correct datatype for HMP_TRANPLANT_MEASURE_%i_WAVELENGTH  --> Default E mean used!", i);
    return eMeanDefault;
  }


  // ef: can remove this?
  if (lambda < 150. || lambda > 230.) {
    LOG(warn) << Form("Wrong value for HMP_TRANPLANT_MEASURE_%i_WAVELENGTH --> Default E mean used!", i);
    return eMeanDefault;
  }

  // find photon energy E in eV from radiation wavelength λ in nm
  nm2eV = 1239.842609;     // 1239.842609 from nm to eV
  photEn = nm2eV / lambda; // photon energy
  return photEn;
}

double HMPIDDCSProcessor::dpVector2Double(const std::vector<DPCOM>& dpVec,
                                          const char* dpString, int i)
{

  double dpVal;
  if (dpVec.size() == 0) {
    LOG(warn) << Form(
      "No Data Point values for HMP_TRANPLANT_MEASURE_%s,%i  "
      "---> Default E mean used!",
      dpString, i);
    return eMeanDefault;
  }

  DPCOM dp = dpVec[0];

  if (dp.id.get_type() == DeliveryType::DPVAL_DOUBLE) {
    dpVal = o2::dcs::getValue<double>(dp);
  } else {
    LOG(warn) << Form(
      "Not correct datatype for HMP_TRANPLANT_MEASURE_%s,%i  "
      "-----> Default E mean used!",
      dpString, i);
    return eMeanDefault;
  }
  return dpVal;
}

bool HMPIDDCSProcessor::evalCorrFactor(double dRefArgon, double dCellArgon,
                                       double dRefFreon, double dCellFreon,
                                       double dPhotEn, int i)
{
  // evaluate correction factor to calculate trasparency (Ref. NIMA 486 (2002)
  // 590-609)

  // Double_t aN1 = AliHMPIDParam::NIdxRad(photEn,tRefCR5);
  // Double_t aN2 = AliHMPIDParam::NMgF2Idx(photEn);
  // Double_t aN3 = 1;                              // Argon Idx

  // Double_t aR1               = ((aN1 - aN2)*(aN1 - aN2))/((aN1 + aN2)*(aN1 +
  // aN2)); Double_t aR2               = ((aN2 - aN3)*(aN2 - aN3))/((aN2 +
  // aN3)*(aN2 + aN3)); Double_t aT1               = (1 - aR1); Double_t aT2 =
  // (1 - aR2); Double_t aCorrFactor       = (aT1*aT1)/(aT2*aT2);

  // evaluate 15 mm of thickness C6F14 Trans

  aConvFactor = 1.0 - 0.3 / 1.8;

  if (dRefFreon * dRefArgon > 0) {
    aTransRad = TMath::Power((dCellFreon / dRefFreon) /
                               (dCellArgon / dRefArgon) * aCorrFactor[i],
                             aConvFactor);
  } else { // ef: removed if mVerbose
    LOGP(warn, "dRefFreon*dRefArgon<0 --> Default E mean used! dRefFreon = {} | dRefArgon = {}", dRefFreon, dRefArgon);    
    return false;
  }

  // evaluate 0.5 mm of thickness SiO2 Trans

  // TMath : Double_t  Exp(Double_t x)
  aTransSiO2 = TMath::Exp(-0.5 / o2::hmpid::Param::lAbsWin(dPhotEn));

  // evaluate 80 cm of thickness Gap (low density CH4) transparency
  aTransGap = TMath::Exp(-80. / o2::hmpid::Param::lAbsGap(dPhotEn));

  // evaluate CsI quantum efficiency
  // if dPhotEn < 6.07267, the value is zero
  aCsIQE = o2::hmpid::Param::qEffCSI(dPhotEn);

  // evaluate total convolution of all material optical properties
  aTotConvolution = aTransRad * aTransSiO2 * aTransGap * aCsIQE;

  sEnergProb += aTotConvolution * dPhotEn;
  double sProbPrev = sProb;
  sProb += aTotConvolution;

  return true;
}
// end Calculate mean photon energy=============================================

// ==== Functions that are called after run is finished ========================

// will return nullptr if there is no entry in Environment-pressure
// DPCOM-vector dpVecEnvPress
std::unique_ptr<TF1> HMPIDDCSProcessor::finalizeEnvPressure()
{
  std::unique_ptr<TF1> pEnv;
  if (dpVecEnv.size() != 0) {

    envPrFirstTime = getMinTime(dpVecEnv);
    envPrLastTime = getMaxTime(dpVecEnv);

    int cntEnvPressure = 0;
    std::unique_ptr<TGraph> pGrPenv;
    pGrPenv.reset(new TGraph);

    for (DPCOM dp : dpVecEnv) {
      auto dpVal = o2::dcs::getValue<double>(dp);
      auto time = dp.data.get_epoch_time(); //-envPrFirstTime
      pGrPenv->SetPoint(cntEnvPressure++, time, dpVal);
    }
    // envPrLastTime -= envPrFirstTime;
    // envPrFirstTime = 0;
 
    // pEnv should not be checked for if nullptr yet, because it will only be defined if the else if or else in the following block is executed:
    if(cntEnvPressure <= 0){
      LOGP(warn, "No entries in Environment Pressure");
    } else if(pGrPenv == nullptr){
      LOGP(warn, "NullPtr in Environment Pressure");
    } else if (cntEnvPressure == 1) {
      pGrPenv->GetPoint(0, xP, yP);
      pEnv.reset(
        new TF1("Penv", Form("%f", yP), envPrFirstTime, envPrLastTime));
    } else {
      // envPrLastTime -= envPrFirstTime;
      // envPrFirstTime = 0;
      pEnv.reset(new TF1("Penv", "1000+x*[0]", envPrFirstTime, envPrLastTime));
      pGrPenv->Fit("Penv", "Q");
    }
    
    // ef: returning nullptr is fine! because,
    // it is checked for in the finalize()
     
    return pEnv;
  }
  LOG(warn) << Form("No entries in environment pressure Penv");
  return pEnv;
}

// returns nullptr if the element in array of DPCOM-vector has no entries
std::unique_ptr<TF1> HMPIDDCSProcessor::finalizeChPressure(int iCh)
{

  std::unique_ptr<TF1> pCh;
  if (dpVecCh[iCh].size() != 0) {
    cntChPressure = 0;

    std::unique_ptr<TGraph> pGrP;
    pGrP.reset(new TGraph);

    chPrFirstTime = getMinTime(dpVecCh[iCh]);
    chPrLastTime = getMaxTime(dpVecCh[iCh]);

    for (DPCOM dp : dpVecCh[iCh]) {
      auto dpVal = o2::dcs::getValue<double>(dp);
      auto time = dp.data.get_epoch_time(); //- chPrFirstTime
      pGrP->SetPoint(cntChPressure++, time, dpVal);
    }
    // chPrLastTime -= chPrFirstTime;
    // chPrFirstTime = 0;
    
    // ef: can we instead of using the smart-pointer arrays, simply call the functions and return the smart-pointer for each iteration?
    

    // ef: would this lead to problems based on the 
    if(pGrP == nullptr || cntChPressure <= 0){
      LOG(warn) << Form("nullptr in chamber-pressure for Pch%i", iCh);
      
    } else if (cntChPressure == 1) {
      pGrP->GetPoint(0, xP, yP);
      (pCh).reset(new TF1(Form("P%i", iCh), Form("%f", yP), chPrFirstTime,
                          chPrLastTime));
      // ef: have to  check for nullptr, because we not can dereference nullptr
      //  (so not equivalent to the environment-pressure)
      if(pCh != nullptr){
        pArrCh[iCh] = *(pCh.get());
      } else{
        LOG(warn) << Form("nullptr in chamber-pressure for Pch%i", iCh);
      }

    } else {
      (pCh).reset(new TF1(Form("P%i", iCh), "[0] + x*[1]", chPrFirstTime,
                          chPrLastTime));
      pGrP->Fit(Form("P%i", iCh), "Q");

      // ef: have to  check for nullptr, because we not can dereference nullptr
      //  (so not equivalent to the environment-pressure)
      if(pCh != nullptr){
        pArrCh[iCh] = *(pCh.get());
      } else{
        LOG(warn) << Form("nullptr in chamber-pressure for P%i", iCh);
      }
    }
    return pCh;
  }
  LOG(warn) << Form("no entries in chamber-pressure for P%i", iCh);

  return pCh;
}

// process Tempout
bool HMPIDDCSProcessor::finalizeTempOut(int iCh, int iRad)
{
  if (dpVecTempOut[3 * iCh + iRad].size() != 0) {
    cntTOut = 0;

    auto minTime = getMinTime(dpVecTempOut[3 * iCh + iRad]);
    auto maxTime = getMaxTime(dpVecTempOut[3 * iCh + iRad]);

    std::unique_ptr<TGraph> pGrTOut;
    pGrTOut.reset(new TGraph);
    
    for (DPCOM dp : dpVecTempOut[3 * iCh + iRad]) {
      auto dpVal = o2::dcs::getValue<double>(dp);
      auto time = dp.data.get_epoch_time(); // -minTime
      pGrTOut->SetPoint(cntTOut++, time, dpVal);
    }
    // maxTime -= minTime;
    // minTime = 0;
    std::unique_ptr<TF1> pTout;
    pTout.reset(
      new TF1(Form("Tout%i%i", iCh, iRad), "[0]+[1]*x", minTime, maxTime));

    // ef: here CCDB-entry will be nullptr if if-case is true
    // this is ok, because on the receiving side of the object, we will check
    // if the entry is nullptr, and in that case use a default value 
    if(pTout == nullptr || pGrTOut == nullptr || cntTOut <= 0 ){
      LOGP(warn, "NullPtr in Temperature out Tout{}{}", iCh, iRad);
      return false;
    } else if (cntTOut == 1) {
      pGrTOut->GetPoint(0, xP, yP);
      pTout->SetParameter(0, yP);
      pTout->SetParameter(1, 0);
    } else {
      pGrTOut->Fit(Form("Tout%i%i", iCh, iRad), "R");
    }

    //ef: in this case everything is ok, and we can set title and put entry in CCDB:
    //    have to check all for nullptr, because pTout is defined before if/else block,
    //    and thus will not be nullptr no matter the outcome of the if/else block 
    if(pTout != nullptr && pGrTOut != nullptr && cntTOut > 0 ){
      pTout->SetTitle(Form("Temp-Out Fit Chamber%i Radiator%i; Time [ms];Temp [C]", iCh, iRad));
      arNmean[6 * iCh + 2 * iRad + 1] = *(pTout.get());
      return true;
    } else{
      LOGP(warn, "NullPtr in Temperature out Tout{}{}", iCh, iRad);
      return false;
    }    
  } else {
    LOGP(warn, "No entries in Temperature out Tout{}{}", iCh, iRad);
    return false;
  }
}



// process Tempin
bool HMPIDDCSProcessor::finalizeTempIn(int iCh, int iRad)
{
  if (dpVecTempIn[3 * iCh + iRad].size() != 0) {
    cntTin = 0;

    auto minTime = getMinTime(dpVecTempIn[3 * iCh + iRad]);
    auto maxTime = getMaxTime(dpVecTempIn[3 * iCh + iRad]);

    std::unique_ptr<TGraph> pGrTIn;
    pGrTIn.reset(new TGraph);

    for (DPCOM dp : dpVecTempIn[3 * iCh + iRad]) {
      auto dpVal = o2::dcs::getValue<double>(dp);
      auto time = dp.data.get_epoch_time(); //-minTime
      pGrTIn->SetPoint(cntTin++, time, dpVal);
    }
    // maxTime -= minTime;
    // minTime = 0;
    std::unique_ptr<TF1> pTin;
    pTin.reset(
      new TF1(Form("Tin%i%i", iCh, iRad), "[0]+[1]*x", minTime, maxTime));

    // ef: here CCDB-entry will be nullptr if if-case is true
    // this is ok, because on the receiving side of the object, we will check
    // if the entry is nullptr, and in that case use a default value 
    if(pTin == nullptr || pGrTIn == nullptr || cntTOut <= 0 ){
      LOGP(warn, "NullPtr in Temperature in Tin{}{}", iCh, iRad);
      return false;
    } else if (cntTOut == 1) {
      pGrTIn->GetPoint(0, xP, yP);
      pTin->SetParameter(0, yP);
      pTin->SetParameter(1, 0);
    } else {
      pGrTIn->Fit(Form("Tin%i%i", iCh, iRad), "R");
    }

    //ef: in this case everything is ok, and we can set title and put entry in CCDB:
    //    have to check all for nullptr, because pTout is defined before if/else block,
    //    and thus will not be nullptr no matter the outcome of the if/else block 
    if(pTin != nullptr && pGrTIn != nullptr && cntTOut > 0 ){
      pTin->SetTitle(Form("Temp-In Fit Chamber%i Radiator%i; Time [ms];Temp [C]", iCh, iRad));
      arNmean[6 * iCh + 2 * iRad + 1] = *(pTin.get());
      return true;
    } else{
      LOGP(warn, "NullPtr in Temperature in Tin{}{}", iCh, iRad);
      return false;
    }    
  } else {
    LOGP(warn, "No entries in Temperature in Tin{}{}", iCh, iRad);
    return false;
  }
}

// returns nullptr if the element in array of DPCOM-vector has no entries
std::unique_ptr<TF1> HMPIDDCSProcessor::finalizeHv(int iCh, int iSec)
{
  std::unique_ptr<TF1> pHvTF;
  if (dpVecHV[3 * iCh + iSec].size() != 0) {

    std::unique_ptr<TGraph> pGrHV;
    pGrHV.reset(new TGraph);
    cntHV = 0;

    hvFirstTime = getMinTime(dpVecHV[3 * iCh + iSec]);
    hvLastTime = getMaxTime(dpVecHV[3 * iCh + iSec]);

    for (DPCOM dp : dpVecHV[3 * iCh + iSec]) {
      auto dpVal = o2::dcs::getValue<double>(dp);
      auto time = dp.data.get_epoch_time(); //- hvFirstTime
      pGrHV->SetPoint(cntHV++, time, dpVal);
    }
    // hvLastTime -= hvFirstTime;
    // hvFirstTime = 0;
    

    // ef: can not check pHvTF for nullptr, because it is defined based on the 
    // the outcome of the following if-block
    if(pGrHV == nullptr || cntHV <= 0){
      LOG(warn) << Form("nullptr in High Voltage for HV%isec%i", iCh, iSec);
    } else if (cntHV == 1) {
      pGrHV->GetPoint(0, xP, yP);


      (pHvTF).reset(new TF1(Form("HV%i_%i", iCh, iSec), Form("%f", yP),
                            hvFirstTime, hvLastTime));
      if(pHvTF!=nullptr){
        pArrHv[3 * iCh + iSec] = *(pHvTF.get());
      } else{
        LOG(warn) << Form("nullptr in High Voltage for HV%i%i", iCh, iSec);
      }

    } else {
      (pHvTF).reset(new TF1(Form("HV%i_%i", iCh, iSec), "[0]+x*[1]",
                            hvFirstTime, hvLastTime));

      if(pHvTF!=nullptr){
        pGrHV->Fit(Form("HV%i_%i", iCh, iSec), "Q");
        pArrHv[3 * iCh + iSec] = *(pHvTF.get());
      } else{
        LOG(warn) << Form("nullptr in High Voltage for HV%i%i", iCh, iSec);
      }
    }
    return pHvTF;
  }
  LOG(warn) << Form("No entries in High Voltage for HV%i%i", iCh, iSec);
  return pHvTF;
}

// Process all arrays of DPCOM-vectors, perform fits
// Called At the end of the Stream of DPs, in
// HMPIDDCSDataProcessorSpec::endOfStream() function
void HMPIDDCSProcessor::finalize()
{
  std::unique_ptr<TF1> pEnv = finalizeEnvPressure();
  for (int iCh = 0; iCh < 7; iCh++) {

    std::unique_ptr<TF1> pChPres = finalizeChPressure(iCh);

    // fills up entries 0..41 of arNmean
    for (int iRad = 0; iRad < 3; iRad++) {
      
      // 6*iCh + 2*iRad 
      bool isTempInValid = finalizeTempIn(iCh, iRad);

      if(isTempInValid == false){
        LOGP(warn, "Tin{}{} not valid! Setting default TF1!", iCh, iRad);
        // this means that the entry was not valid, and thus the vector is not filled
	std::unique_ptr<TF1> pTinDefault;
	pTinDefault.reset(new TF1());
        setDefault(pTinDefault.get(), true);

        // ef: set flag in invalid object, such that it can be read in receiving
 	// side (Ckov reconstruction) as invalid and thus use default value
        arNmean[6 * iCh + 2 * iRad] = *(pTinDefault.get());
      }


      // 6*iCh + 2*iRad + 1
      bool isTempOutValid = finalizeTempOut(iCh, iRad);

      if(isTempOutValid == false){
        LOGP(warn, "Tout{}{} not valid! Setting default TF1!", iCh, iRad);
        // this means that the entry was not valid, and thus the vector is not filled
        std::unique_ptr<TF1> pToutDefault;
	pToutDefault.reset(new TF1());
        setDefault(pToutDefault.get(), true);
        
        // ef: set flag in invalid object, such that it can be read in receiving
 	// side (Ckov reconstruction) as invalid and thus use default value
        arNmean[6 * iCh + 2 * iRad + 1] = *(pToutDefault.get());
      }

    }

    // Fill entries in arQthre
    for (int iSec = 0; iSec < 6; iSec++) {

      std::unique_ptr<TF1> pHV = finalizeHv(iCh, iSec);

      // print the pointers (to see if nullptr)
      LOGP(info, "Finalize ch={} sec={}, pEnv={} pChPres={} pHV={}", iCh, iSec, (void*)pEnv.get(), (void*)pChPres.get(), (void*)pHV.get());

      // only fill if envP, chamP and HV datapoints are all fetched
      if (pEnv != nullptr && pChPres != nullptr && pHV != nullptr) {
        const char* hvChar = (pHV.get())->GetName();     // pArrHv[3 * iCh + iSec]
        const char* chChar = (pChPres.get())->GetName(); //pArrCh[iCh]
        const char* envChar = (pEnv.get())->GetName();
        const char* fFormula =
          "3*10^(3.01e-3*%s - 4.72)+170745848*exp(-(%s+%s)*0.0162012)";
        const char* fName = Form(fFormula, hvChar, chChar, envChar);

        auto minTime = std::max({envPrFirstTime, hvFirstTime, chPrFirstTime});
        auto maxTime = std::min({envPrLastTime, hvLastTime, chPrLastTime});

        std::unique_ptr<TF1> pQthre;
        pQthre.reset(new TF1(Form("HMP_QthreC%iS%i", iCh, iSec), fName, minTime,
                             maxTime));
        pQthre->SetTitle(Form(
          "Charge-Threshold Ch%iSec%i; Time [mS]; Threshold ", iCh, iSec));

        arQthre[6 * iCh + iSec] = *(pQthre.get());

      } else {
        std::unique_ptr<TF1> pQthreDefault;
	pQthreDefault.reset(new TF1());
        setDefault(pQthreDefault.get(), true);
        //setDefault(TF1* f, bool v)// const {f->SetBit(kDefault, v);}
        
        // ef: set flag in invalid object, such that it can be read in receiving
 	// side (Ckov reconstruction) as invalid and thus use default value
        arQthre[6 * iCh + iSec] = *(pQthreDefault.get());
        LOGP(warn, "Missing entry in Charge-Threshold arQthre{} in chamber{} sector{} ", 6 * iCh + iSec, iCh, iSec);

	// print which one is missing:
        if(pEnv == nullptr){
	  LOG(warn) << "Missing entries for environment-pressure Penv";
        }
        if(pChPres == nullptr){
          LOGP(warn, "Missing entries for chamber-pressure P{}", iCh, iCh);
        }
        if(pHV == nullptr){
          LOGP(warn, "Missing entries for High Voltage HV{}{}", iCh, iSec);
        }
      }
    }
  }

  LOG(info) << "======================================== ";

  double eMean = procTrans();

  std::unique_ptr<TF1> pPhotMean;
  pPhotMean.reset(new TF1("HMP_PhotEmean", Form("%f", eMean), 0,
                          mTimeEMean.last - mTimeEMean.first));

  pPhotMean->SetTitle(Form("HMP_PhotEmean; Time [mS]; Photon Energy [eV]"));

  if(pPhotMean != nullptr){
    arNmean[42] = *(pPhotMean.get());
  } else {
    std::unique_ptr<TF1> pPhotMeanDefault;
    setDefault(pPhotMeanDefault.get(), true);
    arNmean[42] = *(pPhotMeanDefault.get());
  }

  // Check entries of CCDB-objects
  checkEntries(arQthre, arNmean);

  // prepare CCDB: =============================================================

  std::map<std::string, std::string> md;
  md["responsible"] = "Erlend Flatland";

  // Refractive index (T_out, T_in, mean photon energy);
  o2::calibration::Utils::prepareCCDBobjectInfo(
    arNmean, mccdbRefInfo, "HMP/Calib/RefIndex", md, mStartValidity, mStartValidity + 3 * o2::ccdb::CcdbObjectInfo::DAY);

  // charge threshold;
  o2::calibration::Utils::prepareCCDBobjectInfo(
    arQthre, mccdbChargeInfo, "HMP/Calib/ChargeCut", md, mStartValidity, mStartValidity + 3 * o2::ccdb::CcdbObjectInfo::DAY);
}

uint64_t HMPIDDCSProcessor::processFlags(const uint64_t flags,
                                         const char* pid)
{

  // function to process the flag. the return code zero means that all is fine.
  // anything else means that there was an issue

  // for now, I don't know how to use the flags, so I do nothing

  if (flags & DPVAL::KEEP_ALIVE_FLAG) {
    LOG(debug) << "KEEP_ALIVE_FLAG active for DP " << pid;
  }
  if (flags & DPVAL::END_FLAG) {
    LOG(debug) << "END_FLAG active for DP " << pid;
  }
  return 0;
}

// extract string from DPID, convert char at specified element in string
// to int
int HMPIDDCSProcessor::aliasStringToInt(const DPID& dpid, std::size_t startIndex)
{

  const std::string inputString(dpid.get_alias());
  char stringPos = inputString[startIndex];
  int charInt = ((int)stringPos) - ((int)'0');
  if (charInt < 10 && charInt >= 0) {
    return charInt;
  } else {
    return -1;
  }
}
} // namespace o2::hmpid
