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
#include <TTree.h>
#include <cassert>

#include "FairLogger.h"
#include "Field/MagneticField.h"
#include "Field/MagFieldFast.h"
#include "TOFBase/Geo.h"

#include "SimulationDataFormat/MCTruthContainer.h"

#include "DetectorsBase/Propagator.h"

#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CommonConstants/GeomConstants.h"
#include "DetectorsBase/GeometryManager.h"

#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include "DataFormatsParameters/GRPObject.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/TrackLTIntegral.h"
#include "ReconstructionDataFormats/TrackHMP.h"

#include "DataFormatsHMP/Digit.h"

#include "GlobalTracking/MatchHMP.h"

#include "TPCBase/ParameterGas.h"
#include "TPCBase/ParameterElectronics.h"
#include "TPCReconstruction/TPCFastTransformHelperO2.h"

#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "HMPIDBase/Param.h"
#include "HMPIDReconstruction/Recon.h"

#include "CommonDataFormat/InteractionRecord.h"

using namespace o2::globaltracking;
using evGIdx = o2::dataformats::EvIndex<int, o2::dataformats::GlobalTrackID>;
using Cluster = o2::hmpid::Cluster;
using Recon = o2::hmpid::Recon;
using MatchInfo = o2::dataformats::MatchInfoHMP;
using MLinfoHMP = o2::dataformats::MLinfoHMP;

using MLinfoHMP = o2::dataformats::MLinfoHMP;
using Trigger = o2::hmpid::Trigger;
using GTrackID = o2::dataformats::GlobalTrackID;
using TrackHMP = o2::dataformats::TrackHMP;
using timeEst = o2::dataformats::TimeStampWithError<float, float>;

ClassImp(MatchHMP);
//==================================================================================================================================================
void MatchHMP::run(const o2::globaltracking::RecoContainer& inp)
{
  ///< running the matching

  mRecoCont = &inp;
  mStartIR = inp.startIR;

  for (int i = 0; i < o2::globaltracking::MatchHMP::trackType::SIZE; i++) {
    mTrackGid[i].clear();
    mMatchedTracks[i].clear();
    mOutHMPLabels[i].clear();
    mTracksWork[i].clear();
    mMatchedTracksIndex[i].clear();
  }

  for (int it = 0; it < o2::globaltracking::MatchHMP::trackType::SIZE; it++) {
    mTracksIndexCache[it].clear();
    if (mMCTruthON) {
      mTracksLblWork[it].clear();
    }
  }

  mHMPTriggersIndexCache.clear();

  bool isPrepareHMPClusters = prepareHMPClusters();

  if (!isPrepareHMPClusters) { // check cluster before of tracks to see also if MC is required
    return;
  }

  // mExtraTPCFwdTime.clear();

  mTimerTot.Start();
  if (!prepareTracks()) {
    return;
  }

  if (mIsITSTPCused || mIsTPCTRDused || mIsITSTPCTRDused || mIsITSTPCTOFused || mIsTPCTOFused || mIsITSTPCTRDTOFused || mIsTPCTRDTOFused) {
    doMatching();
  }

  mIsTPCused = false;
  mIsITSTPCused = false;
  mIsTPCTRDused = false;
  mIsITSTPCTRDused = false;
  mIsTPCTOFused = false;
  mIsITSTPCTRDTOFused = false;
  mIsTPCTRDTOFused = false;
}
//==================================================================================================================================================
bool MatchHMP::prepareTracks()
{

  auto creator = [this](auto& trk, GTrackID gid, float time0, float terr) {
    const int nclustersMin = 0;
    if constexpr (isTPCTrack<decltype(trk)>()) {
      if (trk.getNClusters() < nclustersMin) {
        return true;
      }
      if (std::abs(trk.getQ2Pt()) > mMaxInvPt) {
        return true;
      }
      this->addTPCSeed(trk, gid, time0, terr);
    }
    if constexpr (isTPCITSTrack<decltype(trk)>()) {
      if (trk.getParamOut().getX() < o2::constants::geom::XTPCOuterRef - 1.) {
        return true;
      }
      this->addITSTPCSeed(trk, gid, time0, terr);
    }
    if constexpr (isTRDTrack<decltype(trk)>()) {
      this->addTRDSeed(trk, gid, time0, terr);
    }
    if constexpr (isTPCTOFTrack<decltype(trk)>()) {
      this->addTPCTOFSeed(trk, gid, time0, terr);
    }
    return true;
  };

  mRecoCont->createTracksVariadic(creator);




  for (int it = 0; it < o2::globaltracking::MatchHMP::trackType::SIZE; it++) {
    mMatchedTracksIndex[it].resize(mTracksWork[it].size());
    std::fill(mMatchedTracksIndex[it].begin(), mMatchedTracksIndex[it].end(), -1); // initializing all to -1
  }

  // Unconstrained tracks
  /*
    if (mIsTPCused) {

      // sort tracks in each sector according to their time (increasing in time)
      //  for (int sec = o2::constants::math::NSectors; sec--;) {
      auto& indexCache = mTracksIndexCache[o2::globaltracking::MatchHMP::trackType::UNCONS];
      LOG(debug) << indexCache.size() << " tracks";
      if (!indexCache.size()) {
        return false;
      }
      std::sort(indexCache.begin(), indexCache.end(), [this](int a, int b) {
        auto& trcA = mTracksWork[o2::globaltracking::MatchHMP::trackType::UNCONS][a].second;
        auto& trcB = mTracksWork[o2::globaltracking::MatchHMP::trackType::UNCONS][b].second;
        return ((trcA.getTimeStamp() - trcA.getTimeStampError()) - (trcB.getTimeStamp() - trcB.getTimeStampError()) < 0.);
      });
      // } // loop over tracks of single sector
    } // unconstrained tracks
  */
  // Constrained tracks

  if (mIsITSTPCused || mIsTPCTRDused || mIsITSTPCTRDused || mIsITSTPCTOFused || mIsTPCTOFused || mIsITSTPCTRDTOFused || mIsTPCTRDTOFused) {

    auto& indexCache = mTracksIndexCache[o2::globaltracking::MatchHMP::trackType::CONSTR];
    LOG(debug) << indexCache.size() << " tracks";
    if (!indexCache.size()) {
      return false;
    }
    std::sort(indexCache.begin(), indexCache.end(), [this](int a, int b) {
      auto& trcA = mTracksWork[o2::globaltracking::MatchHMP::trackType::CONSTR][a].second;
      auto& trcB = mTracksWork[o2::globaltracking::MatchHMP::trackType::CONSTR][b].second;
      return ((trcA.getTimeStamp() - mSigmaTimeCut * trcA.getTimeStampError()) - (trcB.getTimeStamp() - mSigmaTimeCut * trcB.getTimeStampError()) < 0.);
    });
    //  } // loop over tracks of single sector
  } // constrained tracks


  LOGP(info, "MatchHMP::prepareTracks() ended ok");
  return true;
}
//______________________________________________
void MatchHMP::addITSTPCSeed(const o2::dataformats::TrackTPCITS& _tr, o2::dataformats::GlobalTrackID srcGID, float time0, float terr)
{
  mIsITSTPCused = true;

  auto trc = _tr.getParamOut();
  o2::track::TrackLTIntegral intLT0 = _tr.getLTIntegralOut();

  timeEst ts(time0, terr);

  addConstrainedSeed(trc, srcGID, ts);
}
//______________________________________________
void MatchHMP::addTRDSeed(const o2::trd::TrackTRD& _tr, o2::dataformats::GlobalTrackID srcGID, float time0, float terr)
{
  if (srcGID.getSource() == o2::dataformats::GlobalTrackID::TPCTRD) {
    mIsTPCTRDused = true;
  } else if (srcGID.getSource() == o2::dataformats::GlobalTrackID::ITSTPCTRD) {
    mIsITSTPCTRDused = true;
  } else if (srcGID.getSource() == o2::dataformats::GlobalTrackID::TPCTRDTOF) {
    mIsTPCTRDTOFused = true;
  } else if (srcGID.getSource() == o2::dataformats::GlobalTrackID::ITSTPCTRDTOF) {
    mIsTPCTRDTOFused = true;
  } else { // shouldn't happen
    LOG(error) << "MatchHMP::addTRDSee: srcGID.getSource() = " << int(srcGID.getSource()) << " not allowed; expected ones are: " << int(o2::dataformats::GlobalTrackID::TPCTRD) << " and " << int(o2::dataformats::GlobalTrackID::ITSTPCTRD) << " and " << int(o2::dataformats::GlobalTrackID::TPCTRDTOF) << " and " << int(o2::dataformats::GlobalTrackID::ITSTPCTRDTOF);
  }

  auto trc = _tr.getOuterParam();

  o2::track::TrackLTIntegral intLT0 = _tr.getLTIntegralOut();

  // o2::dataformats::TimeStampWithError<float, float>
  timeEst ts(time0, terr + mExtraTimeToleranceTRD);

  addConstrainedSeed(trc, srcGID, ts);
}
//______________________________________________
void MatchHMP::addTPCTOFSeed(const o2::dataformats::TrackTPCTOF& _tr, o2::dataformats::GlobalTrackID srcGID, float time0, float terr)
{
  if (srcGID.getSource() == o2::dataformats::GlobalTrackID::TPCTOF) {
    mIsTPCTOFused = true;
  } else if (srcGID.getSource() == o2::dataformats::GlobalTrackID::TPCTRDTOF) {
    mIsTPCTRDTOFused = true;
  } else if (srcGID.getSource() == o2::dataformats::GlobalTrackID::ITSTPCTRDTOF) {
    mIsITSTPCTRDTOFused = true;
  } else { // shouldn't happen
    LOG(error) << "MatchHMP::addTPCTOFCSeed: srcGID.getSource() = " << int(srcGID.getSource()) << " not allowed; expected ones are: " << int(o2::dataformats::GlobalTrackID::TPCTOF) << " and " << int(o2::dataformats::GlobalTrackID::TPCTRDTOF) << " and " << int(o2::dataformats::GlobalTrackID::ITSTPCTRDTOF);
  }

  auto trc = _tr.getParamOut();

  timeEst ts(time0, terr + mExtraTimeToleranceTOF);

  addConstrainedSeed(trc, srcGID, ts);
}
//______________________________________________
void MatchHMP::addConstrainedSeed(o2::track::TrackParCov& trc, o2::dataformats::GlobalTrackID srcGID, timeEst timeMUS)
{
  std::array<float, 3> globalPos;

  // current track index
  int it = mTracksWork[o2::globaltracking::MatchHMP::trackType::CONSTR].size();

  auto prop = o2::base::Propagator::Instance();
  float bxyz[3];
  prop->getFieldXYZ(trc.getXYZGlo(), bxyz);
  double bz = -bxyz[2];

  float pCut = 0.;

  if (TMath::Abs(bz - 5.0) < 0.5) {
    pCut = 0.450;
  }

  if (TMath::Abs(bz - 2.0) < 0.5) {
    pCut = 0.200;
  }

  if (trc.getP() > pCut && TMath::Abs(trc.getTgl()) < 0.544 && TMath::Abs(trc.getPhi() - TMath::Pi()) > (TMath::Pi() * 0.5)) {
    // create working copy of track param
    mTracksWork[o2::globaltracking::MatchHMP::trackType::CONSTR].emplace_back(std::make_pair(trc, timeMUS));

    mTrackGid[o2::globaltracking::MatchHMP::trackType::CONSTR].emplace_back(srcGID);

    if (mMCTruthON) {
      mTracksLblWork[o2::globaltracking::MatchHMP::trackType::CONSTR].emplace_back(mRecoCont->getTPCITSTrackMCLabel(srcGID));
    }

    mTracksIndexCache[o2::globaltracking::MatchHMP::trackType::CONSTR].push_back(it);
  }
}
//______________________________________________
void MatchHMP::addTPCSeed(const o2::tpc::TrackTPC& _tr, o2::dataformats::GlobalTrackID srcGID, float time0, float terr)
{
  mIsTPCused = true;

  std::array<float, 3> globalPos;

  // current track index
  int it = mTracksWork[o2::globaltracking::MatchHMP::trackType::UNCONS].size();

  // create working copy of track param
  timeEst timeInfo;
  // set
  float extraErr = 0;

  auto trc = _tr.getOuterParam();

  float trackTime0 = _tr.getTime0() * mTPCTBinMUS;

  timeInfo.setTimeStampError((_tr.getDeltaTBwd() + 5) * mTPCTBinMUS + extraErr);
  // mExtraTPCFwdTime.push_back((_tr.getDeltaTFwd() + 5) * mTPCTBinMUS + extraErr);

  timeInfo.setTimeStamp(trackTime0);

  trc.getXYZGlo(globalPos);

  mTracksWork[o2::globaltracking::MatchHMP::trackType::UNCONS].emplace_back(std::make_pair(trc, timeInfo));

  if (mMCTruthON) {
    mTracksLblWork[o2::globaltracking::MatchHMP::trackType::UNCONS].emplace_back(mRecoCont->getTPCTrackMCLabel(srcGID));
  }

  mTracksIndexCache[o2::globaltracking::MatchHMP::trackType::UNCONS].push_back(it);
}
//==================================================================================================================================================
bool MatchHMP::prepareHMPClusters()
{
  mHMPClustersArray = mRecoCont->getHMPClusters();
  mHMPTriggersArray = mRecoCont->getHMPClusterTriggers();

  mHMPClusLabels = mRecoCont->getHMPClustersMCLabels();
  mMCTruthON = mHMPClusLabels && mHMPClusLabels->getNElements();

  mNumOfTriggers = 0;

  mHMPTriggersWork.clear();

  int nTriggersInCurrentChunk = mHMPTriggersArray.size();
  LOG(debug) << "nTriggersInCurrentChunk = " << nTriggersInCurrentChunk;
  mNumOfTriggers += nTriggersInCurrentChunk;
  mHMPTriggersWork.reserve(mHMPTriggersWork.size() + mNumOfTriggers);
  for (int it = 0; it < nTriggersInCurrentChunk; it++) {
    const Trigger& clOrig = mHMPTriggersArray[it];
    // create working copy of track param
    mHMPTriggersWork.emplace_back(clOrig);
    //  cache work track index
    mHMPTriggersIndexCache.push_back(mHMPTriggersWork.size() - 1);

    LOG(info) << "mHMPTriggersIndexCache it = " << it;
  }

  // sort hmp events according to their time (increasing in time)
  auto& indexCache = mHMPTriggersIndexCache;
  LOG(debug) << indexCache.size() << " HMP triggers";
  if (!indexCache.size()) {
    return false;
  }

  std::sort(indexCache.begin(), indexCache.end(), [this](int a, int b) {
  auto& clA = mHMPTriggersWork[a];
  auto& clB = mHMPTriggersWork[b];
  auto timeA = o2::InteractionRecord::bc2ns(clA.getBc(), clA.getOrbit());
  auto timeB = o2::InteractionRecord::bc2ns(clB.getBc(), clB.getOrbit());
  return (timeA - timeB) < 0.; });

  return true;
}
//==================================================================================================================================================




// 
void MatchHMP::doMatching()
{


  o2::globaltracking::MatchHMP::trackType type = o2::globaltracking::MatchHMP::trackType::CONSTR;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT; // material correction method


  std::unique_ptr<Recon> recon; 
  recon.reset(new o2::hmpid::Recon());


  //777o2::hmpid::Param* pParam = o2::hmpid::Param::instance();
  // ef moved to init in h
  /*std::unique_ptr<o2::hmpid::Param> pParam;
  pParam->reset(o2::hmpid::Param::instance());*/ 

  o2::hmpid::Param* pParam = o2::hmpid::Param::instance();

  const float kdRadiator = 10.; // distance between radiator and the plane

  //< do the real matching
  auto& cacheTriggerHMP = mHMPTriggersIndexCache; // array of cached HMP triggers indices; reminder: they are ordered in time!
  auto& cacheTrk = mTracksIndexCache[type];       // array of cached tracks indices;
  int nTracks = cacheTrk.size(), nHMPtriggers = cacheTriggerHMP.size();

  LOG(info) << " ************************number of tracks: " << nTracks << ", number of HMP triggers: " << nHMPtriggers;
  if (!nTracks || !nHMPtriggers) {
    return;
  }

  for(const auto& c : cacheTriggerHMP) {
    Printf("const auto& c : cacheTriggerHMP %d", c);
  }
  auto prop = o2::base::Propagator::Instance();

  float bxyz[3];

  double cluLORS[2] = {0};

  LOG(debug) << "Trying to match %d tracks" << cacheTrk.size();

  float timeFromTF = o2::InteractionRecord::bc2ns(mStartIR.bc, mStartIR.orbit);

  
   
  for (int iEvent = 0; iEvent < cacheTriggerHMP.size(); iEvent++) { // events loop

	  auto indexEvent = cacheTriggerHMP[iEvent];
	  Printf("Event number %d   %d ",iEvent, indexEvent);


    auto& event = mHMPTriggersWork[cacheTriggerHMP[iEvent]];
    auto evtTime = o2::InteractionRecord::bc2ns(event.getBc(), event.getOrbit()); // event(trigger) time in ns

	


		const int firstEntry = event.getFirstEntry();
		const int lastEntry = event.getLastEntry();

		LOGP(info, "1st entry last entry {}", mHMPClustersArray.size());
		LOGP(info, "mHMPClustersArray Size  {}", mHMPClustersArray.size());
		for (size_t  j = firstEntry; j <= lastEntry; j++) { // event clusters loop

			if( j >= mHMPClustersArray.size()) {
				
			}

			//const auto& cluster = (o2::hmpid::Cluster&)mHMPClustersArray[j];
			auto& cluster = (o2::hmpid::Cluster&)mHMPClustersArray[j];



			//oneEventClusters.push_back(cluster);
		}

	
    auto& event = mHMPTriggersWork[cacheTriggerHMP[ievt]];
    auto evtTime = event.getIr().differenceInBCMUS(mStartIR);

    int evtTracks = 0;



		// accoutn for 25C -- 20 C in simulation; 
		// dont understand which value to use ..
		// 1.29197 if using mRefIdx(1.28947) + .0005 *(25-5)
		// 1.29044 if using nIdxRad(double eV = 6.675, double temp = 20)
    // double nmean = pParam->meanIdxRad() - 0.0005; // ef TODO: get this from calibration
    double nmean = (1.29197 + 1.29044 )/2.; // ef TODO: get this from calibration
    nmean = 1.2928  - 0.0025; //  ef ; changed ti this to match CkovToolsSingle
    

    
    // collecting all pgoton energies and taking mean : eMean = 6.82 ikke 6.675...
    
	    //nF = 1.2928  - 0.0025; // ef got this from 1 run .. assuming T = 20 for sim 
	     // 1.2928 hvis T = 20; 1.2903 hvis T = 25
	    // np.sqrt(1 + 0.554 * (1239.84 / eV) * (1239.84 / eV) / ((1239.84 / eV) * (1239.84 / eV) - 5769)) - 0.0005 * (temp - 20)
	       
       
    //auto mlEvent = std::make_unique<HmpMLVector>(&oneEventClusters, iEvent, nmean); // ef: initialize as int of Event and clusters relevant to the event

    for (int itrk = 0; itrk < cacheTrk.size(); itrk++) { // tracks loop

    	//oneEventClusters.resize(mHMPClustersArray.size()); this would initialize Cluster objects, dont want to do that
 
      auto& trackWork = mTracksWork[type][cacheTrk[itrk]];
      auto& trackGid = mTrackGid[type][cacheTrk[itrk]];
      auto& trefTrk = trackWork.first;


		  auto indexEvent = cacheTriggerHMP[iEvent];
  	  //Printf("Event %d  Track %d ",iEvent, cacheTrk[itrk]);


      prop->getFieldXYZ(trefTrk.getXYZGlo(), bxyz);
      double bz = -bxyz[2];

      double timeUncert = trackWork.second.getTimeStampError();

      float minTrkTime = (trackWork.second.getTimeStamp() - mSigmaTimeCut * timeUncert);
      float maxTrkTime = (trackWork.second.getTimeStamp() + mSigmaTimeCut * timeUncert);

      if (evtTime < (maxTrkTime + timeFromTF) && evtTime > (minTrkTime + timeFromTF)) {
				Printf("==================================== NEW TRACK ==================");
      // if (evtTime < (maxTrkTime + timeFromTF) && evtTime > (minTrkTime + timeFromTF)) {
      if (evtTime < maxTrkTime && evtTime > minTrkTime) {

        evtTracks++;



        auto hmpTrk = std::make_unique<TrackHMP>(trefTrk); // create a hmpid track to be used for propagation and matching

        //TrackHMP* hmpTrkConstrained = nullptr;    // create a hmpid track to be used for propagation and matching

        std::unique_ptr<TrackHMP> hmpTrkConstrained;

        hmpTrk->set(trefTrk.getX(), trefTrk.getAlpha(), trefTrk.getParams(), trefTrk.getCharge(), trefTrk.getPID());


        // ef: check for simulation if this gets the correct PID::


        double xPc, yPc, xRa, yRa, theta, phi;

        Int_t iCh = intTrkCha(&trefTrk, xPc, yPc, xRa, yRa, theta, phi, bz); // find the intersected chamber for this track
        
        
        //Printf("intTrkCha  xPc %.2f, yPc %.2f, xRa %.2f, yRa %.2f", xPc, yPc, xRa, yRa);


        // LOGP(info, "MatchHMP.cxx Event {}: Track {} Chamber {}" , iEvent, itrk, iCh);
        if (iCh <= 0 || iCh >= 7) {

          continue;
        } // no intersection at all, go next track

        auto matching = std::make_unique<o2::dataformats::MatchInfoHMP>(999999, mTrackGid[type][cacheTrk[itrk]]);
		
        // can these be mveed after if (iCh < 0) {/* statemetns?
        matching->setHMPIDtrk(0, 0, 0, 0);            // no intersection found
        matching->setHMPIDmip(0, 0, 0, 0);            // store mip info in any case
        matching->setIdxHMPClus(99, 99999);           // chamber not found, mip not yet considered
        matching->setHMPsignal(Recon::kNotPerformed); // ring reconstruction not yet performed
        matching->setIdxTrack(trackGid);


        matching->setUnconstrainedPc(xPc, yPc);
        matching->setHMPIDtrk(xRa, yRa, xPc, yPc, theta, phi);
        matching->setIdxHMPClus(iCh, 9999);          // set chamber, index of cluster + cluster size


        MatchInfo matching(999999, mTrackGid[type][cacheTrk[itrk]]);

        matching.setHMPIDtrk(0, 0, 0, 0);            // no intersection found
        matching.setHMPIDmip(0, 0, 0, 0);            // store mip info in any case
        matching.setIdxHMPClus(99, 99999);           // chamber not found, mip not yet considered
        matching.setHMPsignal(Recon::kNotPerformed); // ring reconstruction not yet performed
        matching.setIdxTrack(trackGid);
        TrackHMP hmpTrk(trefTrk); // create a hmpid track to be used for propagation and matching

        hmpTrk.set(trefTrk.getX(), trefTrk.getAlpha(), trefTrk.getParams(), trefTrk.getCharge(), trefTrk.getPID());

        double xPc, yPc, xRa, yRa, theta, phi;

        Int_t iCh = intTrkCha(&trefTrk, xPc, yPc, xRa, yRa, theta, phi, bz); // find the intersected chamber for this track
        if (iCh < 0) {
          continue;
        } // no intersection at all, go next track

        matching.setHMPIDtrk(xPc, yPc, theta, phi); // store initial infos
        matching.setIdxHMPClus(iCh, 9999);          // set chamber, index of cluster + cluster size

        int index = -1;

        double dmin = 999999; //, distCut = 1.;

        bool isOkDcut = kFALSE;
        bool isOkQcut = kFALSE;
        bool isMatched = kFALSE;

        const o2::hmpid::Cluster* bestHmpCluster = nullptr;

        
        auto cluInde = 0;


        LOGP(info, "Track Intersected :::: clusters loop {} -- {}", event.getFirstEntry(), event.getLastEntry());

		    std::vector<Cluster> oneEventClusters;
    		//oneEventClusters.clear();
    		//oneEventClusters.reserve(mHMPClustersArray.size());


				int eventIDClu = -1; // eventID from clusters 
				int eventID = -1; // eventID from clusters 
				LOGP(info, "mHMPClustersArray Size  {}", mHMPClustersArray.size());
				Printf("==================  clusters loop ==================");
        for (int j = event.getFirstEntry(); j <= event.getLastEntry(); j++) { // event clusters loop

	  			if( j >= mHMPClustersArray.size()) {
	    			
          }

          //const auto& cluster = (o2::hmpid::Cluster&)mHMPClustersArray[j];
          auto& cluster = (o2::hmpid::Cluster&)mHMPClustersArray[j];

					eventIDClu = cluster.getEventNumber();

					eventIDClu = cluster.getEventNumber();

<<<<<<< HEAD
          if (cluster.ch() != iCh) {
            // LOGP(info, "cluster.ch() != iCh");
            continue;
          }


	  			int i = j - event.getFirstEntry();

          oneEventClusters.push_back(cluster); //  //this is the important



	  // ef: must loop over oneEventClusters here to make sure these are fulfilled: 
          double qthre = pParam->qCut(); // ef : TODO add chargeCut from calibration!

          if (cluster.q() < 150.) { // 150
          	Printf("Charge too low %.2f", cluster.q());
            continue;
          }


          if (cluster.size() < 3) { // 150
          	Printf("size too low %.2f", cluster.size());
          if (cluster.q() < 150. || cluster.size() > 10) {
            continue;
          }

          // ef: must also check intersection with chambers : oneEventClusters
          // store vector of valid Cluster-indexes per track

          isOkQcut = kTRUE;

          cluLORS[0] = cluster.x();
          cluLORS[1] = cluster.y(); // get the LORS coordinates of the cluster

          double dist = 0.;

          if (TMath::Abs((xPc - cluLORS[0]) * (xPc - cluLORS[0]) + (yPc - cluLORS[1]) * (yPc - cluLORS[1])) > 0.0001) {

            dist = TMath::Sqrt((xPc - cluLORS[0]) * (xPc - cluLORS[0]) + (yPc - cluLORS[1]) * (yPc - cluLORS[1]));
          }

          if (dist < dmin) {
            dmin = dist;
            index =  i;
            // index = oneEventClusters.size() - 1; // not valid w resize/
            bestHmpCluster = &cluster;
            
          }
          LOGP(info, "eventLoop");

        } // event clusters loop


        LOGP(info, "Finished  event clusters loop");
        Printf("==================   ==================");
        // 2. Propagate track to the MIP cluster using the central method

        if (!bestHmpCluster) {/*
          delete hmpTrk;
          hmpTrk = nullptr;
          delete hmpTrkConstrained;
          hmpTrkConstrained = nullptr; */

        if (!bestHmpCluster) {
          oneEventClusters.clear();
          continue;
        }


        matching->setMipClusQ(bestHmpCluster->q());
        
        
        double Dist = TMath::Sqrt((xPc - bestHmpCluster->x()) * (xPc - bestHmpCluster->x()) + (yPc - bestHmpCluster->y()) * (yPc - bestHmpCluster->y()));

=======
        if (!bestHmpCluster) {
          oneEventClusters.clear();
          continue;
        }

>>>>>>> upstream/dev
        TVector3 vG = pParam->lors2Mars(iCh, bestHmpCluster->x(), bestHmpCluster->y());
        float gx = vG.X();
        float gy = vG.Y();
        float gz = vG.Z();
        float alpha = TMath::ATan2(gy, gx);
        float radiusH = TMath::Sqrt(gy * gy + gx * gx);
        if (!(hmpTrk.rotate(alpha))) {
          continue;
        }

        if (!prop->PropagateToXBxByBz(*hmpTrk, radiusH, o2::base::Propagator::MAX_SIN_PHI, o2::base::Propagator::MAX_STEP, matCorr)) {          
	        /*delete hmpTrk;
          hmpTrk = nullptr;
          delete hmpTrkConstrained;
          hmpTrkConstrained = nullptr; */
        if (!prop->PropagateToXBxByBz(hmpTrk, radiusH, o2::base::Propagator::MAX_SIN_PHI, o2::base::Propagator::MAX_STEP, matCorr)) {
          oneEventClusters.clear();
          continue;
        }

        // 3. Update the track with MIP cluster (Improved angular and position resolution - to be used for Cherenkov angle calculation)

        o2::track::TrackParCov trackC(hmpTrk);

        std::array<float, 2> trkPos{0, gz};
        std::array<float, 3> trkCov{0.1 * 0.1, 0., 0.1 * 0.1};

        // auto chi2 = trackC.getPredictedChi2(trkPos, trkCov);
        trackC.update(trkPos, trkCov);

        // 4. Propagate back the constrained track to the radiator radius

        hmpTrkConstrained.reset(new TrackHMP(trackC));
        hmpTrkConstrained->set(trackC.getX(), trackC.getAlpha(), trackC.getParams(), trackC.getCharge(), trackC.getPID());
        if (!prop->PropagateToXBxByBz(*hmpTrkConstrained, radiusH - kdRadiator, o2::base::Propagator::MAX_SIN_PHI, o2::base::Propagator::MAX_STEP, matCorr)) {

	        /*
          delete hmpTrk;
          hmpTrk = nullptr;
          delete hmpTrkConstrained;
          hmpTrkConstrained = nullptr;*/
          continue; 
        TrackHMP hmpTrkConstrained(trackC);
        hmpTrkConstrained.set(trackC.getX(), trackC.getAlpha(), trackC.getParams(), trackC.getCharge(), trackC.getPID());
        if (!prop->PropagateToXBxByBz(hmpTrkConstrained, radiusH - kdRadiator, o2::base::Propagator::MAX_SIN_PHI, o2::base::Propagator::MAX_STEP, matCorr)) {
          oneEventClusters.clear();
          continue;
        }

        float hmpMom = hmpTrkConstrained.getP() * hmpTrkConstrained.getSign();

        matching.setHmpMom(hmpMom);

        // 5. Propagation in the last 10 cm with the fast method


						
        double xPc0 = 0., yPc0 = 0.;
        
        double thetaConst, phiConst;
        intTrkCha(iCh, hmpTrkConstrained.get(), xPc0, yPc0, xRa, yRa, thetaConst, phiConst, bz);
        
        Printf("5. intTrkCha   xRa %.2f xPc0 %.2f yRa %.2f,  yPc0 %.2f",xRa, xPc0, yRa, yPc0);

        // the pc values seem ok until here??
        
        intTrkCha(iCh, &hmpTrkConstrained, xPc0, yPc0, xRa, yRa, theta, phi, bz);

        // 6. Set match information

        int cluSize = bestHmpCluster->size();
        matching->setHMPIDmip(bestHmpCluster->x(), bestHmpCluster->y(), bestHmpCluster->q(), 0); // store mip info in any case
        matching->setMipClusSize(bestHmpCluster->size());
        matching->setIdxHMPClus(iCh, index + 1000 * cluSize); // set chamber, index of cluster + cluster size


        matching->setMipClusPDG(bestHmpCluster->getPDG()); // ef: set event number from cluster

				matching->setMipClusEvent(bestHmpCluster->getEventNumber()); // ef: set event number from cluster


				matching->setMipClusCharge(bestHmpCluster->q()); // ef: set event number from cluster

        //matching->setMIPindex(index); // the position of the MIP in the qrray of clusters



				// ef : wrong to use xPc, yPc here? should be xPc0, yPc0??			
        matching->setHMPIDtrk(xRa, yRa, xPc0, yPc0, thetaConst, phiConst);
        //matching->setHMPIDtrk(xRa, yRa, xPc, yPc, theta, phi);

        

        
        matching->setHMPsignal(pParam->kMipQdcCut);
        matching.setHMPIDmip(bestHmpCluster->x(), bestHmpCluster->y(), (int)bestHmpCluster->q(), 0); // store mip info in any case
        matching.setMipClusSize(bestHmpCluster->size());
        matching.setIdxHMPClus(iCh, index + 1000 * cluSize); // set chamber, index of cluster + cluster size
        matching.setHMPIDtrk(xPc, yPc, theta, phi);

        if (!isOkQcut) {/*
          delete hmpTrk;
          hmpTrk = nullptr;
          delete hmpTrkConstrained;
          hmpTrkConstrained = nullptr; */
          continue;
        if (!isOkQcut) {
          matching.setHMPsignal(pParam->kMipQdcCut);
        }

        // dmin recalculated

        
        dmin = TMath::Sqrt((xPc - bestHmpCluster->x()) * (xPc - bestHmpCluster->x()) + (yPc - bestHmpCluster->y()) * (yPc - bestHmpCluster->y()));

	const auto maxDistAcc = 1.5;
        if (dmin < maxDistAcc) {
          isOkDcut = kTRUE;
        }
        // isOkDcut = kTRUE; // switch OFF cut

        if (!isOkDcut) {
          matching->setHMPsignal(pParam->kMipDistCut); // closest cluster with enough charge is still too far from intersection
          Printf("dmin too high %.2f", dmin);
          matching.setHMPsignal(pParam->kMipDistCut); // closest cluster with enough charge is still too far from intersection
        }

        if (isOkQcut * isOkDcut) {
          isMatched = kTRUE;
        } // MIP-Track matched !!

        Printf("6. : intTrkCha iCh %d || xRa %.2f , xPc0 %.2f ,xMIP %.2f ||  yRa %.2f  yPc0 %.2f  yMIP %.2f \n\n", iCh, xRa, xPc0 ,bestHmpCluster->x(), yRa, yPc0, bestHmpCluster->y());


        // 7. Calculate the Cherenkov angle

        recon->setImpPC(xPc, yPc);                                            // store track impact to PC
        recon->ckovAngle(matching.get(), oneEventClusters, index, nmean, xRa, yRa); // search for Cerenkov angle of this track



        // should this be befoer or after the if (!isMatched) { ?? TODO: ef:
				matching->setRefIndex(nmean);
				matching->setChamber(iCh);





				// set eventID for track from oneEventCluster

				matching->setMipClusEvent(eventIDClu); 


        
				matching->setEventNumber(indexEvent); // 				matching->setEventNumber(iEvent);

				matching->setDistCut(dmin, maxDistAcc); 


<<<<<<< HEAD
        Printf("6. : intTrkCha iCh %d || xRa %.2f , xPc0 %.2f ,xMIP %.2f ||  yRa %.2f  yPc0 %.2f  yMIP %.2f \n\n", iCh, xRa, xPc0 ,bestHmpCluster->x(), yRa, yPc0, bestHmpCluster->y());
=======
        if (!isMatched) {
          mMatchedTracks[type].push_back(matching);
          oneEventClusters.clear();
          mMatchedTracks[type].push_back(*matching);
          continue;
        } // If matched continue...
	

        Printf("bestHmpCluster  x %.2f y %.2f q %d", bestHmpCluster->x(),bestHmpCluster->y(), bestHmpCluster->q());
        
        matching->setMipX(bestHmpCluster->x()); 
        matching->setMipY(bestHmpCluster->y());
        matching->setMipClusCharge(bestHmpCluster->q());
        matching->setMipClusSize(bestHmpCluster->size());

        
        matching->setMipClusEvent(bestHmpCluster->getEventNumber());
				Printf("eventNumber from clu : MIP %d from track : %d", eventIDClu, indexEvent);

			  Printf("Calling match :eventNumber from clu : %d  from track : %d", matching->getMipClusEvent(), matching->getEvent());

                
        Printf("Track::getMIP  x %.2f y %.2f Q %.2f", matching->getMipX(), matching->getMipY(), matching->getMipClusQ());

				// add value indicating matched properly?
				matching->setMatchTrue();
        mMatchedTracks[type].push_back(*matching);
        
        Printf("MATCHED TRUE : CONSTRAINED  xPc0 %.2f, yPc0 %.2f, xRa %.2f, yRa %.2f : BestMIP %.2f %.2f ", xPc0, yPc0, xRa, yRa,bestHmpCluster->x(),bestHmpCluster->y());

        double nmean = pParam->meanIdxRad();

        // 7. Calculate the Cherenkov angle

        recon->setImpPC(xPc, yPc);                                             // store track impact to PC
        recon->ckovAngle(&matching, oneEventClusters, index, nmean, xRa, yRa); // search for Cerenkov angle of this track

        mMatchedTracks[type].push_back(matching);

        oneEventClusters.clear();


        Printf("MATCHED TRUE : UNCONSTRAINED  xPc %.2f, yPc %.2f ", xPc, yPc);


  			double radThick = 1.5, winThick = 0.5, gapThick = 8.0; 
				auto dx = (radThick/2 + winThick + gapThick) * TMath::Cos(phiConst) * TMath::Tan(thetaConst);
				
				auto dy = (radThick/2 + winThick + gapThick) * TMath::Sin(phiConst) * TMath::Tan(thetaConst);
  
  
				double xPcProp = xRa + dx; 
				double yPcProp = yRa + dy;


				Printf("MATCHED TRUE : CONSTRAINED prop xPcProp %.2f, yPcProp %.2f ", xPcProp, yPcProp);

        Printf("Constrained versus unconstrained :  thetaConst %.2f, phiConst %.2f, theta %.2f, phi %.2f \n\n\n\n", thetaConst, phiConst, theta, phi);

      } // if matching in time

      // LOGP(info," MatchHMP.cxx : end track");
    }   // tracks loop


  }     // events loop
  LOGP(info," MatchHMP.cxx : finished all events");

  delete recon;
  recon = nullptr;
}
//==================================================================================================================================================
int MatchHMP::intTrkCha(o2::track::TrackParCov* pTrk, double& xPc, double& yPc, double& xRa, double& yRa, double& theta, double& phi, double bz)
{
  // Static method to find intersection in between given track and HMPID chambers
  // Arguments: pTrk- ESD track; xPc,yPc- track intersection with PC in LORS [cm]
  // Returns: intersected chamber ID or -1
  std::unique_ptr<TrackHMP> hmpTrk;  
  hmpTrk.reset(new TrackHMP(*pTrk));                                        // create a hmpid track to be used for propagation and matching
  for (Int_t i = o2::hmpid::Param::kMinCh; i <= o2::hmpid::Param::kMaxCh; i++) { // chambers loop
    Int_t chInt = intTrkCha(i, hmpTrk.get(), xPc, yPc, xRa, yRa, theta, phi, bz);
    if (chInt >= 0) {
      return chInt;
    }
  } // chambers loop
  return -1; // no intersection with HMPID chambers
} // IntTrkCha()
//==================================================================================================================================================
int MatchHMP::intTrkCha(int ch, o2::dataformats::TrackHMP* pHmpTrk, double& xPc, double& yPc, double& xRa, double& yRa, double& theta, double& phi, double bz)
{
  // Static method to find intersection in between given track and HMPID chambers
  // Arguments: pTrk- HMPID track; xPc,yPc- track intersection with PC in LORS [cm]
  //   Returns: intersected chamber ID or -1


  // ef moved to init in h

  o2::hmpid::Param* pParam = o2::hmpid::Param::instance();



  //o2::hmpid::Param* pParam = o2::hmpid::Param::instance();
  Double_t p1[3], n1[3];
  pParam->norm(ch, n1);
  pParam->point(ch, p1, o2::hmpid::Param::kRad); // point & norm  for middle of radiator plane
  Double_t p2[3], n2[3];
  pParam->norm(ch, n2);
  pParam->point(ch, p2, o2::hmpid::Param::kPc); // point & norm  for entrance to PC plane

  if (pHmpTrk->intersect(p1, n1, bz) == kFALSE) {
    return -1;
  } // try to intersect track with the middle of radiator
  if (pHmpTrk->intersect(p2, n2, bz) == kFALSE) {
    return -1;
  }
  pParam->mars2LorsVec(ch, n1, theta, phi); // track angles at RAD
  pParam->mars2Lors(ch, p1, xRa, yRa);      // TRKxRAD position
  pParam->mars2Lors(ch, p2, xPc, yPc);      // TRKxPC position

  if (pParam->isInside(xPc, yPc, pParam->distCut()) == kTRUE) {
    return ch;
  }          // return intersected chamber
  return -1; // no intersection with HMPID chambers
} // IntTrkCha()
//==================================================================================================================================================
