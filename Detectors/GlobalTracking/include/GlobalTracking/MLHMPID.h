
#ifndef ALICEO2_MLINFOHMP_H
#define ALICEO2_MLINFOHMP_H

// Include the header for the base class
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Trigger.h"

namespace o2
{
namespace globaltracking
{
using MatchInfoHMP = o2::dataformats::MatchInfoHMP;
// Define the derived class
class MLinfoHMP : public MatchInfoHMP
{
  using Cluster = o2::hmpid::Cluster;
  using GTrackID = o2::dataformats::GlobalTrackID;

 private:
  float mxRa, myRa;

  TVector2 mRadImpact;
  float mRefIndex = 1.27;

 public:
  // Default constructor
  MLinfoHMP() = default;

  /*
  MLinfoHMP(int idxHMPClus,
           GTrackID idxTrack,
           float xmip = 0,
           float ymip = 0,
           float xtrk = 0,
           float ytrk = 0,
           float theta = 0,
           float phi = 0,
           float angle = 0,
           float size = 0,
           int idxPhotClus = 0,
           int hmpqn = 0,
           float xRa = 0,
           float yRa = 0,
           float refIndex = 0)

  : MatchInfoHMP(idxHMPClus, idxTrack, xmip, ymip, xtrk, ytrk, theta, phi, angle, size, idxPhotClus, hmpqn),
    mxRa(xRa),
    myRa(yRa),
    mRefIndex(refIndex)  // Initialize the new member variables
  {
      mRadImpact.Set(mxRa, myRa);
  }*/

  MLinfoHMP(const MatchInfoHMP* baseInstance,
            float xRa = 0,
            float yRa = 0,
            float refIndex = -1)
    : MatchInfoHMP(*baseInstance), // Use the copy constructor of the base class
      mxRa(xRa),
      myRa(yRa),
      mRefIndex(refIndex)
  {
    mRadImpact.Set(mxRa, myRa);
  }

  /*
  double getClusters() const
  {
      return mOneEventCluster;
  } */

  void print() const;

  double getRefIndex() const
  {
    return mRefIndex;
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

  TVector2 getRadImpact() const
  {
    return mRadImpact;
  }
  ClassDefNV(MLinfoHMP, 2);
};

class HmpMLVector
{
  using Cluster = o2::hmpid::Cluster;
  using GTrackID = o2::dataformats::GlobalTrackID;

 private:
  std::vector<MLinfoHMP>* mMlTracks; // store fields of tracks

  std::vector<Cluster>* mOneEventCluster;
  float mRefIndex = 1.27;
  int mIevent;

 public:
  // Default constructor
  // HmpMLVector() = default;

  // ef: delete de
  // HmpMLVector() = delete;

  HmpMLVector(std::vector<Cluster>* oneEventcluster = nullptr,
              int iEvent = -1,
              float refIndex = -1)
    : // Use the copy constructor of the base class
      mOneEventCluster(oneEventcluster),
      mRefIndex(refIndex),
      mIevent(iEvent)
  {
  }

  // Add new member functions or override base class functions here
  // For example:
  std::vector<Cluster>* getClusters() const
  {
    return mOneEventCluster;
  }

  double getRefIndex() const
  {
    return mRefIndex;
  }

  void addTrack(std::unique_ptr<MLinfoHMP> hmpTrack)
  {
    mMlTracks->push_back(*hmpTrack);
  }

  void addTrack(const MLinfoHMP& hmpTrack)
  {
    mMlTracks->push_back(hmpTrack);
  }

  std::vector<MLinfoHMP>* getTracks() const
  {
    return mMlTracks;
  }

  const MLinfoHMP& getTrack(int index) const
  {
    return mMlTracks->at(index);
  }
  /*
  MLinfoHMP getTrack(int index) const
  {
    return mMlTracks[index];
  }*/

  ClassDefNV(HmpMLVector, 2);
};

} // namespace globaltracking
} // namespace o2

#endif
