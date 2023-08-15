
struct ShallowDigit {
    uint16_t mQ;
    uint8_t mX;
    uint8_t mY;

    ShallowDigit(uint16_t q, uint8_t x, uint8_t y) : mQ(q), mX(x), mY(y) {}
};

struct ClusterCandidate {
    

    int mCh = 0;
    double mX = 0., mY = 0.;
    double mQ = 0;
    double mChi2 = 0;
    double mXe = 0., mYe = 0.;
    std::vector<ShallowDigit>* mShallowDigits = nullptr;
    std::vector<std::pair<int,int>>* mCandidateStatusVector = nullptr;

    // Constructor based on the order and types you provided
    ClusterCandidate(int ch, double x, double y, double q, double chi2, 
                     double xe, double ye, std::vector<ShallowDigit>* shallowDigits, 
                     std::vector<std::pair<int,int>>* candidateStatusVector) 
        : mCh(ch), mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), 
          mShallowDigits(shallowDigits), mCandidateStatusVector(candidateStatusVector) {}


    //obj.ch, obj.x, obj.y, obj.q, shallowDigits, obj.chi2, obj.xE, obj.yE, candStatus

    void setDigits(const std::vector<ShallowDigit>& shallowDigits) 
    {
        if(!mShallowDigits) {
            mShallowDigits = new std::vector<ShallowDigit>;
        }
        *mShallowDigits = shallowDigits;
    }

    void addCandidateStatus(int iTrack, int hadronCandidateBit) 
    {
        if(!mCandidateStatusVector) {
            mCandidateStatusVector = new std::vector<std::pair<int,int>>;
        }
        mCandidateStatusVector->emplace_back(iTrack, hadronCandidateBit);
    }

    std::vector<std::pair<int,int>>& getCandidateStatus()
    {    
        if(!mCandidateStatusVector) {
            mCandidateStatusVector = new std::vector<std::pair<int,int>>;
        }
        return *mCandidateStatusVector;
    }
};



for(const auto& clusters : clustersVector) // "events loop"
{ 
    // for this event 
    
    // for this event, find intersected chambers: 


    std::sort(MLinfoHMPVector.begin(), MLinfoHMPVector.end(), [](const MLinfoHMP &a, const MLinfoHMP &b) {
        return a.iCh < b.iCh;
    });
    
    std::sort(clusters.begin(), clusters.end(), [](const Cluster &a, const Cluster &b) {
        return a.iCh < b.iCh;
    });


    
    std::vector<MLinfoHMP> sortedTracks[7];
    // Assign MLinfoHMP objects to corresponding vectors based on iCh value
    for (const auto &obj : MLinfoHMPVector) {
        if (obj.iCh >= 0 && obj.iCh <= 6) {
            sortedTracks[obj.iCh].push_back(obj);
        } else {
            std::cerr << "Warning: iCh value out of expected range: " << obj.iCh << std::endl;
        }
    }


    

    // Assuming the range of iCh values is from 0 to 6 (inclusive)
    std::vector<ClusterCandidate> sortedClusters[7];
    // Assign MLinfoHMP objects to corresponding vectors based on iCh value
    for (const auto &obj : clusters) {
        if (obj.iCh >= 0 && obj.iCh <= 6 && sortedTracks[obj.iCh].size() > 0) {

            // make a light copy of digits, just holding the fields charge, x, y
            std::vector<ShallowDigit> shallowDigits;
            shallowDigits.reserve((obj.pDig)->size());
            std::transform(obj.pDig.begin(), obj.pDig.end(), std::back_inserter(shallowDigits),
            [](const Digit* d) {
                return ShallowDigit(d->getQ(), d->getX(), d->getY());
            });


            std::vector<std::pair<int,int>> candStatus = {{0,0}};
            /*
                ClusterCandidate(int ch, double x, double y, double q, double chi2, 
                     double xe, double ye, std::vector<ShallowDigit>* shallowDigits, 
                     std::vector<std::pair<int,int>>* candidateStatusVector) */
            sortedClusters[obj.iCh].emplace_back({obj.ch, obj.x, obj.y, obj.q, obj.chi2, obj.xE, obj.yE, &shallowDigits, &candStatus});
        } else {
            std::cerr << "Warning: iCh value out of expected range: " << obj.iCh << std::endl;
        }
    }
    

    for(int i = 0; i < 7; i++) {
        
        // check if has more than one track
        if(sortedTracks[i].size() < 1) {
            continue;
        }

        auto& clusterPerChamber = sortedClusters[i];
        for(const auto& track : sortedTracks[i]) {

            // pass clusters (and track) by reference, and add its status per track (adding to candStatus vector )
            evaluateClusterTrack(clusterPerChamber, track);
        }


        // now assigned all types of candidates for all clusters in teh given chamber
        // match sortedClusters[i] with sortedTracks[i] --> 
        // 

    }
}


void evaluateClusterTrack(std::vector<ClusterCandidate> clusterPerChamber, MLinfoHMP track);
{


        const auto iEvent = track->getEvent(); // check it corresponds to entry in loop of events?

        const auto momentum = track->getHmpMom()

        const auto nF = track->getRefIndex(); // ef: aon it only gets the mean value ; TODO: in MatchHMP.cxx get calibration value

        const auto nQ = 1.5787;
        const auto nG = 1.0005;

        const auto& ckovHyps = calcCherenkovHyp(momentum, nF); 

        // TODO: this must be changed??


        //double radParams[7] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum,    randomValue.mass};

        float xRad, yRad, xPc, yPc, th, ph;
        track.getHMPIDtrk(xRad, yRad, xPc, yPc, th, ph);

        float xMip, yMip;
        int q, nph;
        // get MIP pos, charge
        track.getHMPIDmip(x, y, q, nph);

        

        const auto& L = 0.5; // 
        double radParams[7] = {xRad, yRad, L, thetaP, phiP, momentum, track.getPID()}; // ef : TODO add PID to MLinfoHMP?
        
        double refIndexes[3] = {nF, nQ, nG};


        // ef: TODO use MIP to get radius and phi in CkovTools: 
        CkovTools ckovTools(radParams, refIndexes, ckovHyps, occupancy, ckovAngle, eventCnt);

        Printf(" Event%d Track%d  : ckovHyps = <%.3f, %.3f> | <%.3f, %.3f> | <%.3f, %.3f>", eventCnt, ckovTools.getMinCkovPion(),ckovTools.getMaxCkovPion(),ckovTools.getMinCkovKaon(), ckovTools.getMaxCkovKaon(),ckovTools.getMinCkovProton(), ckovTools.getMaxCkovProton(), eventCnt); 

        std::array<int, 4> arrayInfo;
        //numBackgroundPhotons, numFoundActualCkov, numActualCkov, numBackgroundLabeledCkov
        std::array<int, 4> arrayInfo;

        // clusterPerChamber by reference, add element to vector {trackNumber, bitHadronStatus}
        ckovTools.segment(clusterPerChamber, arrayInfo, track.getTrackIndex()); // temp --> mapBins
}