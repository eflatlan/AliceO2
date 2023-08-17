/*
struct ShallowDigit {
  int mMotherTrackId;
  int mSourceId;
  int mEventNumber;
  int mTrackId;
  int mParticlePdg;
  uint16_t mQ = 0;
  uint8_t mCh = 0; // 0xFF indicates invalid digit
  uint8_t mPh = 0;
  uint8_t mX = 0;
  uint8_t mY = 0;
  
  uint16_t getQ() const { return mQ; }
  uint8_t getCh() const { return mCh; }
  uint8_t getPh() const { return mPh; }
  uint8_t getX() const { return mX; }
  uint8_t getY() const { return mY; }


  void setEventNumber (int eventNumber) {mEventNumber = eventNumber;}
  int getEventNumber () const  {return mEventNumber;}


  void setMotherId (int motherTrackId) {mMotherTrackId = motherTrackId;}
  int getMotherId () const  {return mMotherTrackId;}


  int getSourceId () const  {return mSourceId;}
    ShallowDigit(uint16_t q, uint8_t x, uint8_t y, Int_t trackId, Int_t particlePdg) : mQ(q), mX(x), mY(y), trackId(mTrackId),particlePdg(mParticlePdg) {}
}; */ 

struct ClusterCandidate {
    

    int mCh = 0;
    double mX = 0., mY = 0.;
    double mQ = 0;
    double mChi2 = 0;
    double mXe = 0., mYe = 0.;
    std::vector<std::pair<int,int>>* mCandidateStatusVector = nullptr;
    std::vector<o2::hmpid::Cluster::Topology> mTopologyVector = nullptr;

    // Constructor based on the order and types you provided
    ClusterCandidate(int ch, double x, double y, double q, double chi2, 
                     double xe, double ye, std::vector<Topology>* topologyVector, 
                     std::vector<std::pair<int,int>>* candidateStatusVector) 
        : mCh(ch), mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), 
          mTopologyVector(topologyVector), mCandidateStatusVector(candidateStatusVector) {}


    //obj.ch, obj.x, obj.y, obj.q, shallowDigits, obj.chi2, obj.xE, obj.yE, candStatus

    void setDigits(const std::vector<Topology>*& topologyVector) 
    {
        if(!mTopologyVector) {
            mTopologyVector = new std::vector<Topology>;
        }
        *mTopologyVector = topologyVector;
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


std::vector<Cluster>* clusterArr = nullptr;
TTree* tCluster = initializeClusterTree(clusterArr);
// clusterArr now initialized correctly

std::vector<o2::dataformats::MatchInfoHMP>* matchArr = nullptr;
TTree* tMatch = initializeMatchTree(matchArr);



std::vector<o2::MCTrack>* mcArr = nullptr;
TTree* tMcTrack = initializeMCTree(matchArr);


// 

int startIndexTrack = 0;
for(int i = 0; i < trigger->size(); i++) //for(const auto& clusters : clustersVector) // "events loop"
{ 

    std::vector<Cluster> oneEventClusters;
    const int firstEntry = pTgr->getFirstEntry();
    const int lastEntry = pTgr->getLastEntry();
    int eventNumber1 = static_cast<o2::hmpid::Cluster>(clusterArr->at(firstEntry));
    int eventNumberLast = static_cast<o2::hmpid::Cluster>(clusterArr->at(lastEntry));
    if(ceventNumberLast != eventNumber1) {Printf("Eventnumber changed??");} // TODO: throw error? ef:


    for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {      
      const auto& clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));

      if(clu.getEventNumber != eventNumber1) {Printf("Eventnumber changed??");}
      else {
          oneEventClusters.push_back(clu);
      }
    }  

    // find entries in tracksOneEvent which corresponds to correct eventNumber
    std::vector<o2::dataformats::MatchInfoHMP>* tracksOneEvent = readMatch(tMatch, matchArr, startIndexTrack);

    // get MC tracks for given event from mc;
    std::vector<o2::MCTrack>* mcTracks = readMC(mcA, tMcTrack, eventNumber1);




    // for this event 


    std::sort(tracksOneEvent.begin(), tracksOneEvent.end(), [](const MLinfoHMP &a, const MLinfoHMP &b) {
        return a.iCh < b.iCh;
    });
    
    std::sort(oneEventClusters.begin(), oneEventClusters.end(), [](const Cluster &a, const Cluster &b) {
        return a.iCh < b.iCh;
    });


    
    std::vector<MLinfoHMP> sortedTracks[7];
    // Assign MLinfoHMP objects to corresponding vectors based on iCh value
    for (const auto &obj : tracksOneEvent) {
        if (obj.iCh >= 0 && obj.iCh <= 6) {
            sortedTracks[obj.iCh].push_back(obj);
        } else {
            std::cerr << "Warning: iCh value out of expected range: " << obj.iCh << std::endl;
        }
    }


    

    // Assuming the range of iCh values is from 0 to 6 (inclusive)
    std::vector<ClusterCandidate> sortedClusters[7];
    // Assign MLinfoHMP objects to corresponding vectors based on iCh value
    for (const auto &obj : oneEventClusters) {
        if (obj.iCh >= 0 && obj.iCh <= 6 && sortedTracks[obj.iCh].size() > 0) {

            // make a light copy of digits, just holding the fields charge, x, y
            /*std::vector<ShallowDigit> shallowDigits;
            shallowDigits.reserve((obj.pDig)->size());
            std::transform(obj.pDig.begin(), obj.pDig.end(), std::back_inserter(shallowDigits),
            [](const Digit* d) {
                return ShallowDigit(d->getQ(), d->getX(), d->getY(), d->getY(), d->getTrackId(), d->getParticlePdg());
            }); */

            const std::vector<o2::hmpid::Cluster::Topology>& topology = obj.getClusterTopology();  // some info about digits associated w cluster

            std::vector<std::pair<int,int>> candStatus = {{0,0}};
            /*
                ClusterCandidate(int ch, double x, double y, double q, double chi2, 
                     double xe, double ye, std::vector<ShallowDigit>* shallowDigits, 
                     std::vector<std::pair<int,int>>* candidateStatusVector) */
            sortedClusters[obj.iCh].emplace_back({obj.ch, obj.x, obj.y, obj.q, obj.chi2, obj.xE, obj.yE, &topology, &candStatus});
        } else {
            std::cerr << "Warning: iCh value out of expected range: " << obj.iCh << std::endl;
        }
    }
    

    for(int i = 0; i < 7; i++) {
        
        // check if has more than one track --> this means there is no candidates
        if(sortedTracks[i].size() < 1) {
            continue;
        }

        auto& clusterPerChamber = sortedClusters[i];


        std::vector<float> mipCharges;
        // fill charges of MIPs
        for(const auto& track : sortedTracks[i]) {
            float xMip, yMip;
            track.getHMPIDmip(xMip, yMip, q, nph);
            mipCharges.emplace_back(q);
        }

        for(const auto& track : sortedTracks[i]) {

            // pass clusters (and track) by reference, and add its status per track (adding to candStatus vector )

            // for each clusterPerChamber we will have a "candidate-status" for each of the photons, this is a vector of length of sortedTracks[i].size();
            // and holds the fields 

            
            // get MCTrack correspondign to trackId
            const auto mcTrackIndex = track.getTrackIndex();
            const o2::MCTrack& mcTrack = getMCEntry(mcTracks, mcTrackIndex);

            const int mcTrackPdg = mcTrack.GetPdgCode();
            // add mcTrackPdg
            evaluateClusterTrack(clusterPerChamber, track, mipCharges, mcTrackPdg);
        }

        // save ClusterPerChamber

        // save trackInfo : 
        //                  trackIndex : mIdxTrack
        //                  scalar : p, refIndex,
        //                  trackInfo theta, phi, xRad, yRad; {xPc yPc is redundant}
        //                  MIP            : x, y, q
        //                  MCTruth : pdg code of track || pdg code of MIP?
          



        // now assigned all types of candidates for all clusters in teh given chamber
        // match sortedClusters[i] with sortedTracks[i] --> 
        // 

    }
}


void evaluateClusterTrack(std::vector<ClusterCandidate>& clusterPerChamber, const MLinfoHMP& track, const std::vector<float>& mipCharges, int mcTrackPdg);
{

        const auto iEvent = track.getEvent(); // check it corresponds to entry in loop of events?

        const auto momentum = track.getHmpMom();

        const auto nF = track.getRefIndex(); // ef: aon it only gets the mean value ; TODO: in MatchHMP.cxx get calibration value

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

        float xMip, yMip;
        track.getHMPIDmip(xMip, yMip, q, nph);

        

        const auto& L = 0.5; // 
        double radParams[7] = {xRad, yRad, L, thetaP, phiP, momentum, track.getPID()}; // ef : TODO add PID to MLinfoHMP?
        
        double refIndexes[3] = {nF, nQ, nG};


        // ef: TODO use MIP to get radius and phi in CkovTools: 
        CkovTools ckovTools(radParams, refIndexes, ckovHyps, occupancy, ckovAngle, eventCnt);

        Printf(" Event%d Track%d  : ckovHyps = <%.3f, %.3f> | <%.3f, %.3f> | <%.3f, %.3f>", eventCnt, ckovTools.getMinCkovPion(),ckovTools.getMaxCkovPion(),ckovTools.getMinCkovKaon(), ckovTools.getMaxCkovKaon(),ckovTools.getMinCkovProton(), ckovTools.getMaxCkovProton(), eventCnt); 

        std::array<int, 4> arrayInfo;
        //numBackgroundPhotons, numFoundActualCkov, numActualCkov, numBackgroundLabeledCkov
        std::array<int, 4> arrayInfo;



        // add charge of MIPS here? to skip them in candidates ? 

        // clusterPerChamber by reference, add element to vector {trackNumber, bitHadronStatus}


        // mcTrackPdg check that it matches with clusterPDG?
        ckovTools.segment(clusterPerChamber, arrayInfo, track.getTrackIndex(), mipCharges, mipX, mipY, mcTrackPdg); // temp --> mapBins
}