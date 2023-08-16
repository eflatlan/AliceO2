

template <typename T>
void initFileIn(std::unique_ptr<TFile>& mFile, std::unique_ptr<TTree>& mTree, const std::string& filename, const std::string& firstTree, const std::string& secondTree, const std::string& firstBranch, const std::string& secondBranch, , std::vector<T>*& dataFromFile) {
 
  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  // Create the TFIle
  mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());

  mTree.reset((TTree *)mFile->Get(firstTree));
  if (!mTree) {
    mTree.reset((TTree *)mFile->Get(secondTree));
  }

  if (!mTree) {
    LOGP(warn, "HMPID DigitToClusterSpec::initFileIn() : Did not find {} tree in file {}", secondTree, filename.c_str());
    return;
    std::exit(0);
  }

  if ((mTree->GetBranchStatus(firstBranch)) == 1) {
    mTree->SetBranchAddress(firstBranch, &dataFromFile);
  } else if ((mTree->GetBranchStatus(secondBranch)) == 1) {
    mTree->SetBranchAddress(secondBranch, &dataFromFile);
  } else {
   LOG(warn)
        << "HMPID DigitToClusterSpec::initFileIn() : Error in branches!" << endl;
    return;
    std::exit(0);
  }

  //mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);
  mTree->Print("toponly");
}

void readTreeEntries() {
    // Open the ROOT file
    auto matchFile = std::make_unique<TFile>("o2match_hmp.root");
    auto clusterFile = std::make_unique<TFile>("hmpidclus.root");
    auto mcFile = std::make_unique<TFile>("o2sim_Kine.root");

    std::unique_ptr<TTree> matchTree, clusTree, kineTRee; ///< input tree

    std::vector<Clusters>* mClustersFromFile;
    std::vector<o2::dataformats::MLinfoHMP>* mTracksFromFile;
    std::vector<o2::MCTrack>* mCarloFromFile;

    initFileIn(matchFile, matchTree,  "",  "", "HMPMatchInfo", "HMPMatchInfo", &mTracksFromFile );
    /*if (mUseMC) {
        mTree->SetBranchAddress("MatchHMPMCTruth", &mLabelHMPPtr);
    }*/ 

    initFileIn(clusterFile, clusTree,  "o2hmp",  "o2sim", "HMPIDClusters", "HMPIDclusters", &mClustersFromFile);

    initFileIn(mcFile, kineTRee,  "o2hmp",  "o2sim", "MCTrack", "MCTrack", &mCarloFromFile);

    

    /*
    auto currEntry = mTree->GetReadEntry() + 1;
    assert(currEntry < mTree->GetEntries()); // this should not happen
    mTree->GetEntry(currEntry);*/



    // Loop over tree entries
    Long64_t nMatchEvents = matchTree->GetEntries();
    Long64_t nClusterEvents = clusTree->GetEntries();
    Long64_t nKineEvents = kineTRee->GetEntries();
    LOGP(info, "nMatchEvents {}, nClusterEvents {nClusterEvents},  nKineEvents {nKineEvents}");

    // loop over nKineTree see if we have HMP tracks and clusters for the given event

    for (Long64_t iKine = 0; iKine < nKineEvents; iKine++) {
        if (nKineEvents->GetReadEntry() + 1 >= nKineEvents->GetEntries()) {
        }  
        else {
            kineTRee->GetEntry(iKine);  // This fills the above variables with data

            for (int i = 0; i < mCarloFromFile->size(); ++i) {
              const auto& mcTrack = (*mcAmCarloFromFilerr)[i];

              // 
              if (i == trackID) {
                Printf("Particle %d: pdg = %d, pT = %f, px = %f, py = %f, pz = %f, vx = %f, vy = %f, vz = %f", i, mcTrack.GetPdgCode(), TMath::Abs(mcTrack.GetStartVertexMomentumX() * mcTrack.GetStartVertexMomentumX() + mcTrack.GetStartVertexMomentumY() * mcTrack.GetStartVertexMomentumY()), mcTrack.GetStartVertexMomentumX(), mcTrack.GetStartVertexMomentumY(), mcTrack.GetStartVertexMomentumZ(), mcTrack.GetStartVertexCoordinatesX(), mcTrack.GetStartVertexCoordinatesY(), mcTrack.GetStartVertexCoordinatesZ());
              }
            }

            matchTree->GetEntry(iKine);  // This fills the above variables with data
            clusTree->GetEntry(iKine);  // This fills the above variables with data

            // Now you can use the data
            //std::cout << "Entry " << i << ": branchA = " << branchAData << ", branchB = " << branchBData << std::endl;
        }
    }

    //file->Close();




  tKine->GetEntry(eventID);
  for (int i = 0; i < mCarloFromFile->size(); ++i) {
    const auto& mcTrack = (*mCarloFromFile)[i];
    if (i == trackID) {
      Printf("Particle %d: pdg = %d, pT = %f, px = %f, py = %f, pz = %f, vx = %f, vy = %f, vz = %f", i, mcTrack.GetPdgCode(), TMath::Abs(mcTrack.GetStartVertexMomentumX() * mcTrack.GetStartVertexMomentumX() + mcTrack.GetStartVertexMomentumY() * mcTrack.GetStartVertexMomentumY()), mcTrack.GetStartVertexMomentumX(), mcTrack.GetStartVertexMomentumY(), mcTrack.GetStartVertexMomentumZ(), mcTrack.GetStartVertexCoordinatesX(), mcTrack.GetStartVertexCoordinatesY(), mcTrack.GetStartVertexCoordinatesZ());
    }
  }

}