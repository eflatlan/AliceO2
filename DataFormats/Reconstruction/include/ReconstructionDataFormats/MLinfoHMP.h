// Include the header for the base class
#include "MatchInfoHMP.h"

namespace o2
{
namespace dataformats
{

// Define the derived class
class MLinfoHMP : public MatchInfoHMP
{ 

private:

std::vector<cluster>* mOneEventCluster;
float mxRa, myRa;
float mRefIndex = 1.27;
TVector2 mRadImpact;

public:
    // Default constructor
    MLinfoHMP() = default;
    
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
                         float RefIndex = 0)

    : MatchInfoHMP(idxHMPClus, idxTrack, xmip, ymip, xtrk, ytrk, theta, phi, angle, size, idxPhotClus, hmpqn),
      mxRa(xRa), 
      myRa(yRa), 
      mRefIndex(RefIndex)  // Initialize the new member variables
    {
        mRadImpact.Set(mxRa, myRa)
    }

    MLinfoHMP(const MatchInfoHMP* baseInstance,                         
                        float xRa = 0, 
                        float yRa = 0, 
                        float RefIndex = 0)
    : MatchInfoHMP(*baseInstance), // Use the copy constructor of the base class 
    mXra(xRa), 
    mYra(yRa), 
    mRefIndex(RefIndex)
    {
        mRadImpact.Set(mxRa, myRa);
    }

    /*
    double getClusters() const 
    {
        return mOneEventCluster;
    } */ 

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

class HmpMLVector
{ 

  private:

    std::vector<MLinfoHMP>* mMlTracks; // store fields of tracks 

    std::vector<cluster>* mOneEventCluster;
    float mxRa, myRa;
    float mRefIndex = 1.27;
    TVector2 mRadImpact;

  public:
      // Default constructor
      HmpMLVector() = default;
      


      HmpMLVector(std::vector<Cluster>* oneEventcluster = {}, 
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

      void addTrack(const MLinfoHMP& hmpTrack) 
      {
        mMlTracks->push_back(hmpTrack)
      }

      
      std::vector<MLinfoHMP>* getTracks() const 
      {
        return mMlTracks;
      }

      MLinfoHMP getTracks(int index) const 
      {
        return mMlTracks[index];
      }

};


} // namespace dataformats
} // namespace o2
