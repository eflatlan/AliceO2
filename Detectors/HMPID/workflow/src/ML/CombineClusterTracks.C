

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

    std::unique_ptr<TTree> matchTree, clusTree; ///< input tree



    std::vector<Clusters>* mClustersFromFile;
    std::vector<o2::dataformats::MLinfoHMP>* mTracksFromFile;

    initFileIn(matchFile, matchTree,  "",  "", "HMPMatchInfo", "HMPMatchInfo", &mTracksFromFile );
    /*if (mUseMC) {
        mTree->SetBranchAddress("MatchHMPMCTruth", &mLabelHMPPtr);
    }*/ 

    initFileIn(clusterFile, clusTree,  "o2hmp",  "o2sim", "HMPIDClusters", "HMPIDclusters", &mClustersFromFile);

    

    /*
    auto currEntry = mTree->GetReadEntry() + 1;
    assert(currEntry < mTree->GetEntries()); // this should not happen
    mTree->GetEntry(currEntry);*/



    // Loop over tree entries
    Long64_t nEntries = matchTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {

        else {
            matchTree->GetEntry(i);  // This fills the above variables with data
            clusTree->GetEntry(i);  // This fills the above variables with data

            // Now you can use the data
            std::cout << "Entry " << i << ": branchA = " << branchAData << ", branchB = " << branchBData << std::endl;
        }
    }

    //file->Close();
}