// Include the header for the base class
#include "MatchInfoHMP.h"

namespace o2
{
namespace dataformats
{

// Define the derived class
class ExtendedMatchInfoHMP : public MatchInfoHMP
{

 private:
  std::vector<cluster>* mOneEventCluster;
  float mxRa, myRa;
  float mRefIndex = 1.27;
  TVector2 mRadImpact;

 public:
  // Default constructor
  ExtendedMatchInfoHMP() = default;

  ExtendedMatchInfoHMP(std::vector<Cluster>* oneEventcluster = {},
                       int iEvent,
                       float refIndex)
    : MatchInfoHMP(*baseInstance), // Use the copy constructor of the base class
      mOneEventCluster(oneEventcluster),
      mRefIndex(refIndex),
      mIevent(iEvent)
  {
    mRadImpact.Set(mxRa, myRa);
  }

  // Add new member functions or override base class functions here
  // For example:
  double getClusters() const
  {
    return mOneEventCluster;
  }

  double getRefIndex() const
  {
    return mRefIndex;
  }

  void getHMPIDtrk(float& xRa, float& yRa, float& xPc, float& yPc, float& th, float& ph) const
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
};

} // namespace dataformats
} // namespace o2
