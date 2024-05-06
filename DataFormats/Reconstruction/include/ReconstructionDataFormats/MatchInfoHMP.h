w// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MatchInfoHMP.h
/// \brief Class to store the output of the matching to HMPID

#ifndef ALICEO2_MATCHINFOHMP_H
#define ALICEO2_MATCHINFOHMP_H

#include "ReconstructionDataFormats/TrackLTIntegral.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "CommonDataFormat/EvIndex.h"
#include <TVector2.h>

namespace o2
{
namespace dataformats
{
class MatchInfoHMP
{
  using GTrackID = o2::dataformats::GlobalTrackID;

 public:
  MatchInfoHMP(int idxHMPClus, GTrackID idxTrack, float xmip = 0, float ymip = 0, float xtrk = 0, float ytrk = 0, float theta = 0, float phi = 0, float angle = 0, float size = 0, int idxPhotClus = 0, int hmpqn = 0, float hmpMom = 0) : mIdxHMPClus(idxHMPClus), mIdxTrack(idxTrack), mCkovAngle(angle), mMipX(xmip), mMipY(ymip), mTrkX(xtrk), mTrkY(ytrk), mTrkTheta(theta), mTrkPhi(phi), mMipCluSize(size), mIdxPhotClus(idxPhotClus), mHMPqn(hmpqn), mHmpMom(hmpMom)
  { // Initialize the mPhotCharge array
    for (int i = 0; i < 10; i++) {
      mPhotCharge[i] = 0.0;
    }
  };
  MatchInfoHMP() = default;

  void setIdxHMPClus(int ch, int idx) { mIdxHMPClus = ch * 1000000 + idx; }
  int getIdxHMPClus() const { return mIdxHMPClus; }

  void setIdxTrack(GTrackID index) { mIdxTrack = index; }
  int getTrackIndex() const { return mIdxTrack.getIndex(); }

  void setMipclusIndex(int index) { mipIndex = index; } // ef : remove these three
  int getMipclusIndex() const { return mipIndex; }
  int mipIndex = -1;

  GTrackID getTrackRef() const { return mIdxTrack; }

  void setMipX(float x) { mMipX = x; }
  float getMipX() const { return mMipX; }

  void setMipY(float y) { mMipY = y; }
  float getMipY() const { return mMipY; }

  void setTrkX(float x) { mTrkX = x; }
  float getTrkX() const { return mTrkX; }

  void setMipClusQ(float q) { mMipCluQ = q; }
  float getMipClusQ() const { return mMipCluQ; }

  void setTrkY(float y) { mTrkY = y; }
  float getTrkY() const { return mTrkY; }

  void setHMPsignal(float angle) { mCkovAngle = angle; }
  
  float getHMPsignal() const
  {
    if (mCkovAngle > 0) {
      return mCkovAngle - (Int_t)mCkovAngle;
    } else {
      return mCkovAngle;
    }
  }


  // ef added
  void setHMPsignalMassHyp(float angle) { mCkovAngleMassHyp = angle; }
  
  float getHMPsignalMassHyp() const
  {
    if (mCkovAngleMassHyp > 0) {
      return mCkovAngleMassHyp - (Int_t)mCkovAngle;
    } else {
      return mCkovAngleMassHyp;
    }
  }

  void setMassHypNumPhot(int nPhot) {
    mNPhotsMassHyp = nPhot;
  }


  int getMassHypNumPhot() const { return mNPhotsMassHyp; }

  //


  void setDistCut(float dist, float distThre)
  {
    mDist = dist;
    mDistThre = distThre;
  }
  void getDistCut(float& dist, float& distThre) const
  {
    dist = mDist;
    distThre = mDistThre;
  }

  void setEventNumberFromTrack(int eventNumberTrack) { mEventNumberTrack = eventNumberTrack; }
  int getEventNumberFromTrack() const { return mEventNumberTrack; }  

  void setMipClusPDG(int pdg) { mMipCluPDG = pdg; }
  int getMipClusEventPDG() const { return mMipCluPDG; }

  void setMipClusEvent(int event) { mMipCluEvent = event; }
  int getMipClusEvent() const { return mMipCluEvent; }

  void setMipClusCharge(float charge) { mMipCluCharge = charge; }
  int getMipClusCharge() const { return mMipCluCharge; }

  void setMipClusSize(int size) { mMipCluSize = size; }
  int getMipClusSize() const { return mMipCluSize; }

  void setNPhots(int n) { mNPhots = n; }
  int getNPhots() const { return mNPhots; }

  void setPhotIndex(int idx) { mIdxPhotClus = idx; }
  int getPhotIndex() const { return mIdxPhotClus; }

  float getOccupancy() const { return (Int_t)mCkovAngle / 10.0; }

  void setHMPIDtrk(float x, float y, float th, float ph)
  {
    mTrkX = x;
    mTrkY = y;
    mTrkTheta = th;
    mTrkPhi = ph;
  }

  void setHMPIDtrk(float xRa, float yRa, float x, float y, float th, float ph)
  {

    mxRa = xRa;
    myRa = yRa;

    mTrkX = x;
    mTrkY = y;
    mTrkTheta = th;
    mTrkPhi = ph;
  }

  int getEvent() const
  {
    return mEvent;
  }

  void getHMPIDrad(float& xRad, float& yRad) const
  {
    xRad = mxRa;
    yRad = myRa;
  }

  void setHMPIDrad(float xRad, float yRad)
  {
    mxRa = xRad;
    myRa = yRad;
  }


  void getHMPIDtrk(float& xRad, float& yRad, float& xPc, float& yPc, float& th, float& ph) const
  {
    xRad = mxRa;
    yRad = myRa;
    xPc = mTrkX;
    yPc = mTrkY;
    th = mTrkTheta;
    ph = mTrkPhi;
  }

  void getHMPIDtrk(float& xPc, float& yPc, float& th, float& ph) const
  {

    xPc = mTrkX;
    yPc = mTrkY;
    th = mTrkTheta;
    ph = mTrkPhi;
  }

  void setRefIndex(float refIndex)
  {
    mRefIndex = refIndex;
  }

  double getRefIndex() const
  {
    return mRefIndex;
  }

  void setChamber(int iChamber)
  {
    miCh = iChamber;
  }

  int getChamber() const
  {
    return miCh;
  }

  void setEventNumber(int iEvent)
  {
    mEvent = iEvent;
  }

  void setHMPIDmip(float x, float y, int q, int nph = 0)
  {
    mMipX = x;
    mMipY = y;
    mHMPqn = 1000000 * nph + q;
  }

  void getHMPIDmip(float& x, float& y, int& q, int& nph) const
  {
    x = mMipX;
    y = mMipY;
    q = mHMPqn % 1000000;
    nph = mHMPqn / 1000000;
  }

  void setParticlePdg(int particlePdg) { mParticlePdg = particlePdg; }
  int getParticlePdg() const { return mParticlePdg; }

  void setHmpMom(float p) { mHmpMom = p; }
  float getHmpMom() const { return mHmpMom; }

  void setPhotCharge(float* chargeArray)
  {
    for (int i = 0; i < 10; i++) {
      mPhotCharge[i] = chargeArray[i];
    }
  }

  const float* getPhotCharge() const { return mPhotCharge; }

  void print() const;

  void setMatchTrue() { isMatched = true; }

  bool getMatchStatus() const { return isMatched; }

 private:
  bool isMatched = false;
  float mxRa, myRa;

  float mRefIndex = 1.27;
  int mEvent;
  int miCh;
  
  int mEventNumberTrack = -1;

  float mDist = -1., mDistThre = -1.; // distance to MIP, cut for distance used

  /*



   void setMipClusPDG(int pdg) { mMipCluPDG = pdg; }
   int getMipClusEventPDG() const { return mMipCluPDG; }

   void setMipClusEvent(int event) { mMipCluEvent = event; }
   int getMipClusEvent() const { return mMipCluEvent; }

   void setMipClusCharge(int size) { mMipCluCharge = size; }
   int getMipClusCharge() const { return mMipCluCharge; }

   void setMipClusSize(int size) { mMipCluSize = size; }
   int getMipClusSize() const { return mMipCluSize; }

 */

 protected:
  float mMipCluQ = 0.0; // MIP cluster charge

  int mMipCluPDG;     // for sim: the PDG code of the MIP matched w track
  int mMipCluEvent;   // for sim: the Evnet of the  MIP matched w track
  int mMipCluCharge;  // for sim: the charge of the MIP matched w track
  int mParticlePdg;   // for sim: the PDG code of the track
  int mIdxHMPClus;    // Idx for HMP cluster
  GTrackID mIdxTrack; // Idx for track
  float mMipX;        // local x coordinate of macthed cluster
  float mMipY;        // local y coordinate of macthed cluster
  float mTrkX;        // local x coordinate of extrapolated track intersection point
  float mTrkY;        // local y coordinate of extrapolated track intersection point
  float mTrkTheta;    // theta track
  float mTrkPhi;      // phi track
  float mCkovAngle;   // emission angle value
  int mMipCluSize;    // MIP cluster size
  int mNPhots;        // number of candidate photo-electrons
  int mIdxPhotClus;   // index of the first photo
  int mHMPqn;         // 1000000*number of photon clusters + QDC
  float mHmpMom;      // track momentum at HMPID chambers
  float mPhotCharge[10] = {};

  float mCkovAngleMassHyp; // ef aadded
  int mNPhotsMassHyp;


  ClassDefNV(MatchInfoHMP, 2);
};
} // namespace dataformats
} // namespace o2
#endif

