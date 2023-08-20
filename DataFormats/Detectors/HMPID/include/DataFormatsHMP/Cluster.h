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

struct Topology {
	uint8_t posX = 0;
	uint8_t posY = 0;
	uint16_t q = 0;
	int pdg = 0;
	int tid = 0;
	int mid = 0;
	int eid = 0;
	int currentCluVecSize = 0;

	Topology() = default;  // Default constructor, uses default member initializers.

  Topology(uint8_t x, uint8_t y, uint16_t q_val, int p, int t, int m, int e, int c)
        : posX(x), posY(y), q(q_val), pdg(p), tid(t), mid(m), eid(e), currentCluVecSize(c) {}
};


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
  /*void getTopologyVector(std::vector<int>& outDigsX, 
                       std::vector<int>& outDigsY, 
                       std::vector<int>& outDigsQ,  
                       std::vector<int>& outDigsPDG, 
                       std::vector<int>& outDigsTID, 
                       std::vector<int>& outDigsMID) 
	{
		  /* outDigsX = digsX;
		  outDigsY = digsY;
		  outDigsQ = digsQ;
		  outDigsPDG = digsPDG;
		  outDigsTID = digsTID;
		  outDigsMID = digsMID;* /
	}*/



  /*
  void setClusterTopology(const std::vector<Topology>& topVector)
  { 



    // <o2::hmpid::Digit*>* mDigs
    for(const auto& dig : topVector){
      const auto& posX = dig.posX; // pos of digit
      const auto& posY = dig.posY ;
      const auto& q = dig.q;
      const auto& pdg = dig.pdg;
      const auto& tid = dig.tid;
      const auto& mid = dig.mid;
      mTopologyVector.emplace_back(Topology{posX, posY, q, pdg, tid, mid});
    }  
    
	}*/ 



  

  void setLastTopologyIndex(int index) { mLastTopIndex = index;} 
  void setFirstTopologyIndex(int index) { mFirstTopIndex = index;} 

  int getLastTopologyIndex() const { return mLastTopIndex;} 
  int getFirstTopologyIndex() const { return mFirstTopIndex;} 

	void setClusterTopology(std::vector<o2::hmpid::Topology>& topVector, const int currCluVecSize)
	{ 
      
		  //std::vector<o2::hmpid::Topology> topVector;
		  // Check for null or empty mDigs
		  if(!mDigs || mDigs->empty()) {
		      setEventNumber(-1);
		      return;// topVector;
		  }

		  // Assuming all digits have the same event number, so take it from the first digit
		  setEventNumber((*mDigs)[0]->getEventNumber());
			
		  // Reserve space for vectors (Optional, but can improve efficiency)
		  const size_t digSize = mDigs->size();

			//topVector.reserve(digSize);

		  for(const auto& dig : *mDigs) {
					if(dig) {			
						const auto& posX = dig->getX(); // pos of digit
						const auto& posY = dig->getY() ;
						const auto& q = dig->getQ();
						const auto& pdg = dig->mParticlePdg;
						const auto& tid = dig->mTrackId;
						const auto& mid = dig->getMotherId();
						const auto& eid = dig->getEventNumber();

						int pdgCat = 0;
						if(pdg == 50000050 || pdg == 50000051) {
							pdgCat = 4;
						} else if (pdg == 211 || pdg == 111) { // Pion
							pdgCat = 3;
						} else if (pdg == 321) { // Kaon
							pdgCat = 2;
						} else if (pdg == 2212) { // Proton
							pdgCat = 1;
						} else  { // Proton
							pdgCat = 0;
						}

      			topVector.emplace_back(o2::hmpid::Topology{posX, posY, q, pdgCat, tid, mid, eid, currCluVecSize});
					}
		  }
			//return topVector;
	}


  Cluster() : mCh(-1), mSi(-1), mSt(kEmp), mBox(-1), mNlocMax(-1), mMaxQpad(-1), mMaxQ(-1), mQRaw(0), mQ(0), mErrQ(-1), mXX(0), mErrX(-1), mYY(0), mErrY(-1), mChi2(-1) {}

  // Methods
  // void draw(Option_t *opt=""); //overloaded TObject::Print() to draw cluster in current canvas
  void print(Option_t* opt = "") const;                                                  // overloaded TObject::Print() to print cluster info
  static void fitFunc(int& iNpars, double* deriv, double& chi2, double* par, int iflag); // fit function to be used by MINUIT
  void cleanPointers()
  {
    mDigs = nullptr;
  }

  void coG();                                                                            // calculates center of gravity
  void corrSin();                                                                        // sinoidal correction
  void digAdd(const o2::hmpid::Digit* pDig);                                             // add new digit to the cluster
  const o2::hmpid::Digit* dig(int i) const { return mDigs ? (*mDigs)[i] : nullptr; }     // pointer to i-th digi
  const std::vector<const o2::hmpid::Digit*>* getDigits() const { return mDigs; }
  void setDigits(std::vector<const o2::hmpid::Digit*>* v = nullptr) { mDigs = v; }
  inline bool isInPc();                                                                  // check if is in the current PC
  void reset();                                                                          // cleans the cluster
  // void setClusterParams(float xL, float yL, int iCh); //Set AliCluster3D part
  int solve(std::vector<o2::hmpid::Cluster>* pCluLst,  std::vector<o2::hmpid::Topology>& pTopVector, float* pSigmaCut, bool isUnfold); // solve cluster: MINUIT fit or CoG
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

  //int  digsXArr[9], digsYArr[9], digsQArr[9], digsPDGArr[9];//, digsTIDArr[9], digsMIDArr[9];//.push_back(dig

  std::vector<int>  digsX, digsY, digsQ, digsPDG;//, digsTID, digsMID;//.push_back(dig->getX()); // pos of digit

	int mLastTopIndex = -1, mFirstTopIndex = -1;
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
  std::vector<const o2::hmpid::Digit*>* mDigs = nullptr; //! list of digits forming this cluster
  //std::vector<Topology> mTopologyVector;

  public  : 
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


  int mEventNumber; // ef: 

 public:
  static bool fgDoCorrSin; // flag to switch on/off correction for Sinusoidal to cluster reco
  void setEventNumber(int eNum) { mEventNumber = eNum; }
  int getEventNumber() const { return mEventNumber; }
  ClassDefNV(Cluster, 4);
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
