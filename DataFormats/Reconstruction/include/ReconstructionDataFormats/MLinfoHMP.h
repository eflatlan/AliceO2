
#ifndef ALICEO2_MLINFOHMP_H
#define ALICEO2_MLINFOHMP_H

// Include the header for the base class
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include <TVector2.h>

namespace o2
{
namespace dataformats
{
using MatchInfoHMP = o2::dataformats::MatchInfoHMP;
// Define the derived class
class MLinfoHMP : public MatchInfoHMP
{

  using GTrackID = o2::dataformats::GlobalTrackID;

 private:
  float mxRa, myRa;

  TVector2 mRadImpact;
  float mRefIndex = 1.27;
  int mEvent;
  int miCh;

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
            int iCh = -1,
            float xRa = -1.0,
            float yRa = -1.0,
            float refIndex = -1.0,
            int event = -1)
    : MatchInfoHMP(*baseInstance), // Use the copy constructor of the base class
      miCh(iCh),
      mxRa(xRa),
      myRa(yRa),
      mRefIndex(refIndex),
      mEvent(event)
  {
    mRadImpact.Set(mxRa, myRa);
  }

  void print() const;

  int getChamber() const
  {
    return miCh;
  }

  int getEvent() const
  {
    return mEvent;
  }

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

} // namespace dataformats
} // namespace o2

#endif
