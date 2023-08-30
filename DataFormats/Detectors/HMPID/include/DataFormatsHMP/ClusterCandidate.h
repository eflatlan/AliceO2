#ifndef CLUSTER_CANDIDATE_H
#define CLUSTER_CANDIDATE_H
#include <TObject.h>
namespace o2
{
namespace hmpid
{
class ClusterCandidate : public TObject { 


  public:
   
    int mSize = 0;
    int mCh = 0;
    double mX = 0., mY = 0.;
    int mQ = 0;
    double mChi2 = 0;
    double mXe = 0., mYe = 0.;
    int mPDG = -1;

    // vector e.l. som holder truth? // i.e., for hver track, set MIP og trackIndex fra track
    int trackId = -1;
    bool isMip = false;

    int mCandidateStatus = 0;// = {{0,0}}; do not initialize
    int trackNumber = 0;// = {{0,0}}; do not initialize
    
    
    
    
    // std::vector<o2::hmpid::Cluster::Topology> mTopologyVector = nullptr;


    ClusterCandidate() = default;
    // ClusterCandidate() : mCh(0), mX(0), mY(0), mQ(0), mChi2(0), mXe(0), mYe(0), mPDG(-1) {}


    // Constructor based on the order and types you provided
    ClusterCandidate(int ch, double x, double y, int q, double chi2, 
                     double xe, double ye, /*std::vector<Topology>* topologyVector,*/ int pdg, int size) 
        : mCh(ch), mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), 
          /*mTopologyVector(topologyVector),*/ mPDG(pdg), mSize(size) {}


    //obj.ch, obj.x, obj.y, obj.q, shallowDigits, obj.chi2, obj.xE, obj.yE, candStatus


    /*
    void setDigits(const std::vector<Topology>*& topologyVector) 
    {
        if(!mTopologyVector) {
            mTopologyVector = new std::vector<Topology>;
        }
        *mTopologyVector = topologyVector;
    } */

    void setCandidateStatus(int hadronCandidateBit) /*const*/
    {

      mCandidateStatus = hadronCandidateBit;

    }

    int getCandidateStatus()
    {   
        return mCandidateStatus;
    }
    ClassDef(ClusterCandidate, 1);

};
}
}
#endif // CLUSTER_CANDIDATE_H

