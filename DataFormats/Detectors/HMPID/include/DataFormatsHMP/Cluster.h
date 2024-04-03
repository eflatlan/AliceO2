// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef DETECTORS_HMPID_BASE_INCLUDE_HMPIDDATAFORMAT_CLUSTER_H_

#define DETECTORS_HMPID_BASE_INCLUDE_HMPIDDATAFORMAT_CLUSTER_H_

#include "CommonDataFormat/InteractionRecord.h"

#include "CommonDataFormat/RangeReference.h"

#include "DataFormatsHMP/Digit.h"

#include "HMPIDBase/Param.h"

#include "FairLogger.h"

namespace o2

{

namespace hmpid

{

/// \class Cluster

/// \brief HMPID cluster implementation

class Cluster

{

 public:
  enum EClusterStatus { kFrm,

                        kCoG,

                        kLo1,

                        kUnf,

                        kMax,

                        kNot,

                        kEdg,

                        kSi1,

                        kNoLoc,

                        kAbn,

                        kBig,

                        kEmp = -1 }; // status flags

 public:
   /*float getEnergy() const { return mEnergy; }

   void setEnergy(float energy) { mEnergy = energy; }

   int getPDG() const { return mParticlePdg; }

   void setPDG(int pdg) { mParticlePdg = pdg; }

   void setEventNumberFromTrack(int eventNumberTrack) { mEventNumberTrack =
   eventNumberTrack; }

   int getEventNumberFromTrack() const { return mEventNumberTrack; }

   void setEventNumber(int eventNumber) { mEventNumber = eventNumber; }

   int getEventNumber() const { return mEventNumber; }

   void setTrackId(int tid) { mTrackId = tid; }

   int getTrackId() const { return mTrackId; }

   void setMotherId(int motherTrackId) { mMotherTrackId = motherTrackId; }

   int getMotherId() const { return mMotherTrackId; }

   int getSourceId() const { return mSourceId; }


   void setInfoFromDigits(const int currCluVecSize)

   {

     //  Check for null or empty mDigs

     if (!mDigs || mDigs->empty()) {

       setEventNumber(-1);

       return;
     }

     std::sort(mDigs->begin(), mDigs->end(), [](const o2::hmpid::Digit* a, const
   o2::hmpid::Digit* b) { return a->getCharge() < b->getCharge();
     });

     const int clux = this->x();

     const int cluy = this->y();

     std::vector<const o2::hmpid::Digit*> digs = *mDigs;

     std::sort(digs.begin(), digs.end(), [clux, cluy](const o2::hmpid::Digit* a,
   const o2::hmpid::Digit* b) { int xa, ya;

       int cha;

       o2::hmpid::Digit::pad2Absolute(a->getPadID(), &cha, &xa, &ya);

       int xb, yb, chb;

       o2::hmpid::Digit::pad2Absolute(b->getPadID(), &chb, &xb, &yb);

       // calc dist to Clu pos

       auto digaR = (xa - clux) * (xa - clux) + (ya - cluy) * (ya - cluy);

       auto digbR = (xb - clux) * (xb - clux) + (yb - cluy) * (yb - cluy);

       return digaR < digbR;
     });

     //setEventNumber((*mDigs)[0]->getEventNumber());

     // assuming main digit now is the first...

     setTrackId((*mDigs)[0]->getTrackId());

     setMotherId((*mDigs)[0]->getMotherId());

     setPDG((*mDigs)[0]->getPDG());


     LOGP(info, "======================");

     LOGP(info, "Based on charge : pdg {} mother {} tid {}",
   (*mDigs)[0]->getPDG(), (*mDigs)[0]->getMotherId(),
   (*mDigs)[0]->getTrackId());

     LOGP(info, "Based on pos : pdg {} mother {} tid {}", digs[0]->getPDG(),
   digs[0]->getMotherId(), digs[0]->getTrackId());



     const size_t digSize = mDigs->size();
   }

   void setInfoFromDigits2(const int currCluVecSize, std::vector<int>
   indicesResolvedDig)

   {

     //  Check for null or empty mDigs

     if (!mDigs || mDigs->empty()) {

       setEventNumber(-1);

       return;
     }

     std::vector<const o2::hmpid::Digit*> digs;
     std::vector<const o2::hmpid::Digit*> digs2;

     for (int digIndex : indicesResolvedDig) {
       if (digIndex > mDigs->size())
         LOGP(error, "digIndex {} mDigs->size() {}", digIndex, mDigs->size());
       else {
         digs.push_back((*mDigs)[digIndex]);
         digs2.push_back((*mDigs)[digIndex]);
       }
     }

     std::sort(digs2.begin(), digs2.end(), [](const o2::hmpid::Digit* a, const
   o2::hmpid::Digit* b) { return a->getCharge() < b->getCharge();
     });

     const int clux = this->x();

     const int cluy = this->y();

     std::sort(digs.begin(), digs.end(), [clux, cluy](const o2::hmpid::Digit* a,
   const o2::hmpid::Digit* b) { int xa, ya;

       int cha;

       o2::hmpid::Digit::pad2Absolute(a->getPadID(), &cha, &xa, &ya);

       int xb, yb, chb;

       o2::hmpid::Digit::pad2Absolute(b->getPadID(), &chb, &xb, &yb);

       // calc dist to Clu pos

       auto digaR = (xa - clux) * (xa - clux) + (ya - cluy) * (ya - cluy);

       auto digbR = (xb - clux) * (xb - clux) + (yb - cluy) * (yb - cluy);

       return digaR < digbR;
     });

     if (digs2.size() == 0) {
       return;
     }
     setEventNumber((digs2)[0]->getEventNumber());

     // assuming main digit now is the first...

     setTrackId((digs2)[0]->getTrackId());

     setMotherId((digs2)[0]->getMotherId());

     setPDG((digs2)[0]->getPDG());

     LOGP(info, "======================");

     LOGP(info, "Based on charge : pdg {} mother {} tid {}",
   (digs2)[0]->getPDG(), (digs2)[0]->getMotherId(), (digs2)[0]->getTrackId());

     // LOGP(info, "Based on pos : pdg {} mother {} tid {}", digs[0]->getPDG(),
   digs[0]->getMotherId(), digs[0]->getTrackId());

     const size_t digSize = mDigs->size();
   }
   */
   Cluster()
       : mCh(-1), mSi(-1), mSt(kEmp), mBox(-1), mNlocMax(-1), mMaxQpad(-1),
         mMaxQ(-1), mQRaw(0), mQ(0), mErrQ(-1), mXX(0), mErrX(-1), mYY(0),
         mErrY(-1), mChi2(-1) {}

   // Methods

   // void draw(Option_t *opt=""); //overloaded TObject::Print() to draw cluster
   // in current canvas

   void print(Option_t *opt = "")
       const; // overloaded TObject::Print() to print cluster info

   static void fitFunc(int &iNpars, double *deriv, double &chi2, double *par,
                       int iflag); // fit function to be used by MINUIT

   void cleanPointers()

   {
     mDigs = nullptr;
   }

  void coG(); // calculates center of gravity

  void corrSin(); // sinoidal correction

  void digAdd(const o2::hmpid::Digit* pDig); // add new digit to the cluster

  const o2::hmpid::Digit* dig(int i) const { return mDigs ? (*mDigs)[i] : nullptr; } // pointer to i-th digi

  const std::vector<const o2::hmpid::Digit*>* getDigits() const { return mDigs; }

  void setDigits(std::vector<const o2::hmpid::Digit*>* v = nullptr) { mDigs = v; }

  inline bool isInPc(); // check if is in the current PC

  void reset(); // cleans the cluster

  // void setClusterParams(float xL, float yL, int iCh); //Set AliCluster3D part

  int solve(std::vector<o2::hmpid::Cluster>* pCluLst, float* pSigmaCut, bool isUnfold); // solve cluster: MINUIT fit or CoG

  // ef :added
  int solveMC(std::vector<o2::hmpid::Cluster>* pCluLst, float* pSigmaCut, bool isTryUnfold, std::map<int, std::vector<int>>& resolvedIndicesMap);
  void findClusterSizeMC(int i, float* pSigmaCut, std::vector<int>& indicesResolved);

  // Getters

  // Setters

  void doCorrSin(bool doCorrSin) { fgDoCorrSin = doCorrSin; } // Set sinoidal correction

  void setX(float x) { mXX = x; }

  void setY(float y) { mYY = y; }

  void setQ(float q)

  {

    mQ = q;

    if (mQ > 4095) {

      mQ = 4095;
    }
  }

  void setQRaw(float qRaw)

  {

    mQRaw = qRaw;

    if (mQRaw > 4095) {

      mQRaw = 4095;
    }
  }

  void setSize(int size) { mSi = size; }

  void setCh(int chamber) { mCh = chamber; }

  void setChi2(float chi2) { mChi2 = chi2; }

  void setStatus(int status) { mSt = status; }

  void findClusterSize(int i, float* pSigmaCut); // Find the clusterSize of deconvoluted clusters

  // public:

protected:
  // int  digsXArr[9], digsYArr[9], digsQArr[9], digsPDGArr[9];//, digsTIDArr[9], digsMIDArr[9];//.push_back(dig

  // std::vector<int>  digsX, digsY, digsQ, digsPDG;//, digsTID, digsMID;//.push_back(dig->getX()); // pos of digit

  int mCh; // chamber number

  int mSi; // size of the formed cluster from which this cluster deduced

  int mSt; // flag to mark the quality of the cluster

  int mBox; // box contaning this cluster

  int mNlocMax; // number of local maxima in formed cluster

  int mMaxQpad; // abs pad number of a pad with the highest charge

  double mMaxQ; // that max charge value

  double mQRaw; // QDC value of the raw cluster

  double mQ; // QDC value of the actual cluster

  double mErrQ; // error on Q

  double mXX; // local x postion, [cm]

  double mErrX; // error on x postion, [cm]

  double mYY; // local y postion, [cm]

  double mErrY; // error on y postion, [cm]

  double mChi2; // some estimator of the fit quality

  std::vector<const o2::hmpid::Digit*>* mDigs = nullptr; //! list of digits forming this cluster

 public:
   // ef : added const methods
   int box() const { return mBox; } // Dimension of the cluster

   int ch() const { return mCh; } // chamber number

   int size() const { return mSi; } // returns number of pads in formed cluster

   int status() const { return mSt; } // Status of cluster

   float qRaw() const { return mQRaw; } // raw cluster charge in QDC channels

   float q() const { return mQ; } // given cluster charge in QDC channels

   float qe() const { return mErrQ; } // Error in cluster charge in QDC channels

   float x() const { return mXX; } // cluster x position in LRS

   float xe() const { return mErrX; } // cluster charge in QDC channels

   float y() const { return mYY; } // cluster y position in LRS

   float ye() const { return mErrY; } // cluster charge in QDC channels

   float chi2() const { return mChi2; } // chi2 of the fit

   int box() { return mBox; } // Dimension of the cluster

   int ch() { return mCh; } // chamber number

   int size() { return mSi; } // returns number of pads in formed cluster

   int status() { return mSt; } // Status of cluster

   float qRaw() { return mQRaw; } // raw cluster charge in QDC channels

   float q() { return mQ; } // given cluster charge in QDC channels

   float qe() { return mErrQ; } // Error in cluster charge in QDC channels

   float x() { return mXX; } // cluster x position in LRS

   float xe() { return mErrX; } // cluster charge in QDC channels

   float y() { return mYY; } // cluster y position in LRS

   float ye() { return mErrY; } // cluster charge in QDC channels

   float chi2() { return mChi2; } // chi2 of the fit

 public:
  static bool fgDoCorrSin; // flag to switch on/off correction for Sinusoidal to cluster reco

  ClassDefNV(Cluster, 5);
};

} // namespace hmpid

namespace framework

{

template <typename T>

struct is_messageable;

template <>

struct is_messageable<o2::hmpid::Cluster> : std::true_type {
};

} // namespace framework

} // namespace o2

#endif /* DETECTORS_HMPID_BASE_INCLUDE_HMPIDDATAFORMAT_CLUSTER_H_ */
