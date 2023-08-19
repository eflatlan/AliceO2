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
struct Topology {
    int diffX;
    int diffY;
    int q;
    int pdg;
    int tid;
    int mid;
};


 public:
  Cluster(const Cluster& other)
    : mCh(other.mCh), mSi(other.mSi), mSt(other.mSt),
      mBox(other.mBox), mNlocMax(other.mNlocMax), mMaxQpad(other.mMaxQpad),
      mMaxQ(other.mMaxQ), mQRaw(other.mQRaw), mQ(other.mQ),
      mErrQ(other.mErrQ), mXX(other.mXX), mErrX(other.mErrX),
      mYY(other.mYY), mErrY(other.mErrY), mChi2(other.mChi2),
      mEventNumber(other.mEventNumber)
  {
      // Copy the mDigs vector, you should deep copy if needed. Here it's shallow copy.
      mDigs = other.mDigs;

      // Copy the mTopologyVector
      mTopologyVector = other.mTopologyVector;
  }

  Cluster() : mCh(-1), mSi(-1), mSt(kEmp), mBox(-1), mNlocMax(-1), mMaxQpad(-1), mMaxQ(-1), mQRaw(0), mQ(0), mErrQ(-1), mXX(0), mErrX(-1), mYY(0), mErrY(-1), mChi2(-1) {}

  //~Cluster() {  }

  // Methods
  // void draw(Option_t *opt=""); //overloaded TObject::Print() to draw cluster in current canvas
  void print(Option_t* opt = "") const;                                                  // overloaded TObject::Print() to print cluster info
  static void fitFunc(int& iNpars, double* deriv, double& chi2, double* par, int iflag); // fit function to be used by MINUIT
  
  void cleanPointers()
  {
		///delete mDigs; 
    //mDigs = nullptr;
  }
  void coG();                                                                            // calculates center of gravity
  void corrSin();                                                                        // sinoidal correction
  void digAdd(const o2::hmpid::Digit& pDig);                                             // add new digit to the cluster
  const o2::hmpid::Digit dig(int i) const { if(mDigs.size() < i) { return mDigs[i];}
  o2::hmpid::Digit dig; return dig;  }     // pointer to i-th digi

  void setClusterTopology()
  { 
    mTopologyVector.reserve(mDigs.size());


    int eventNumber = -1;
    
    if(!mDigs.empty()) {
      eventNumber = (mDigs[0]).getEventNumber();

    }
    setEventNumber(eventNumber);


    // <o2::hmpid::Digit*>* mDigs
    for(const auto& dig : mDigs){
      const auto& posX = dig.getX(); // pos of digit
      const auto& posY = dig.getY();
      const auto& q = dig.getQ();
      const auto& pdg = dig.mParticlePdg;
      const auto& tid = dig.mTrackId;
      const auto& mid = dig.mMotherTrackId;
      mTopologyVector.emplace_back(Topology{posX, posY, q, pdg, tid, mid});
    }
    
  }


  const int getNumDigits() const { return mDigs.size(); }
  const std::vector<o2::hmpid::Digit> getDigits() const { return mDigs; }
  void setDigits(std::vector<o2::hmpid::Digit> v) { mDigs = v; }
  inline bool isInPc();                                                                  // check if is in the current PC
  void reset();                                                                          // cleans the cluster
  // void setClusterParams(float xL, float yL, int iCh); //Set AliCluster3D part
  

  //int solve(std::vector<o2::hmpid::Cluster>* pCluLst, float* pSigmaCut, bool isUnfold); // solve cluster: MINUIT fit or CoG
int solve(/*std::vector<o2::hmpid::Cluster>& pCluLst, */Cluster& clu, float* pSigmaCut, bool isTryUnfold);
int solve(std::vector<o2::hmpid::Cluster>& cluLst,  float* pSigmaCut, bool isTryUnfold);
  // Getters

	/*ef : TODO: why are these not marked const??*/
  int box() { return mBox; }     // Dimension of the cluster
  int ch() { return mCh; }       // chamber number
  int size() { return mSi; }     // returns number of pads in formed cluster
  int status() { return mSt; }   // Status of cluster
  float qRaw() { return mQRaw; } // raw cluster charge in QDC channels
  float q() { return mQ; }       // given cluster charge in QDC channels

  float qe() { return mErrQ; }   // Error in cluster charge in QDC channels
  float x() { return mXX; }      // cluster x position in LRS
  float xe() { return mErrX; }   // cluster charge in QDC channels
  float y() { return mYY; }      // cluster y position in LRS
  float ye() { return mErrY; }   // cluster charge in QDC channels
  float chi2() { return mChi2; } // chi2 of the fit

	int box() const { return mBox; }     // Dimension of the cluster
	int ch() const { return mCh; }       // chamber number
	int size() const { return mSi; }     // returns number of pads in formed cluster
	int status() const { return mSt; }   // Status of cluster
	float qRaw() const { return mQRaw; } // raw cluster charge in QDC channels
	float q() const { return mQ; }       // given cluster charge in QDC channels
	float qe() const { return mErrQ; }   // Error in cluster charge in QDC channels
	float x() const { return mXX; }      // cluster x position in LRS
	float xe() const { return mErrX; }   // cluster charge in QDC channels
	float y() const { return mYY; }      // cluster y position in LRS
	float ye() const { return mErrY; }   // cluster charge in QDC channels
	float chi2() const { return mChi2; } // chi2 of the fit


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

  void setEventNumber(int eNum) { mEventNumber = eNum; }
  int getEventNumber() const { return mEventNumber; }

  void setSize(int size) { mSi = size; }
  void setCh(int chamber) { mCh = chamber; }
  void setChi2(float chi2) { mChi2 = chi2; }
  void setStatus(int status) { mSt = status; }
  void findClusterSize(int i, float* pSigmaCut); // Find the clusterSize of deconvoluted clusters

  const std::vector<Topology>& getTopologyVector() const { return mTopologyVector;}

  
  // public:
 protected:

  std::vector<Topology> mTopologyVector;

//mTopologyVector.emplace_back(diffX, diffY, q, pdg, tid, mid);
  int mCh;                                    // chamber number
  int mSi;                                    // size of the formed cluster from which this cluster deduced
  int mSt;                                    // flag to mark the quality of the cluster
  int mBox;                                   // box contaning this cluster
  int mNlocMax;                               // number of local maxima in formed cluster
  int mMaxQpad;                               // abs pad number of a pad with the highest charge
  double mMaxQ;                               // that max charge value
  double mQRaw;                               // QDC value of the raw cluster
  double mQ;                                  // QDC value of the actual cluster
  double mErrQ;                               // error on Q
  double mXX;                                 // local x postion, [cm]
  double mErrX;                               // error on x postion, [cm]
  double mYY;                                 // local y postion, [cm]
  double mErrY;                               // error on y postion, [cm]
  double mChi2;                               // some estimator of the fit quality
  std::vector<o2::hmpid::Digit> mDigs; //! list of digits forming this cluster

  int mEventNumber; // ef: 

 public:
  static bool fgDoCorrSin; // flag to switch on/off correction for Sinusoidal to cluster reco

  ClassDefNV(Cluster, 3);
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
