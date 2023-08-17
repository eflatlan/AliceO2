#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)

#include "DataFormatsHMP/Cluster.h"
#include "HMPIDReconstruction/Clusterer.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TApplication.h>
#include <TF1.h>
#include <TH2F.h>
#include <TLine.h>
#include <TList.h>
#include <TROOT.h> // gRoot
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>

#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <vector>
#include <fairlogger/Logger.h>


#include "CommonDataFormat/InteractionRecord.h"

// C++ header files and libraries
#include <math.h>
#include <chrono>
#include <thread>
#include <ctime>
#include <fstream>
#include <iostream>
#include <gsl/gsl>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <filesystem>
#include <iostream>
#include <string>
#endif


using std::this_thread::sleep_for;
using std::vector, std::cout, std::cin, std::endl;
using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger, o2::hmpid::Clusterer;


using Clusters = o2::hmpid::Cluster;
using Cluster = o2::hmpid::Cluster;//, o2::hmpid::Digit, o2::hmpid::Trigger, o2::hmpid::Clusterer;

//template <typename T>

void initFileIn(std::unique_ptr<TFile>& mFile, std::unique_ptr<TTree>& mTree, const std::string& firstTree, const std::string& secondTree, const std::string& firstBranch, const std::string& secondBranch) {
 
  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  // Create the TFIle
  //mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());

  mTree.reset((TTree *)mFile->Get(firstTree.c_str()));
  if (!mTree) {
    mTree.reset((TTree *)mFile->Get(secondTree.c_str()));
  }

  if (!mTree) {
    LOGP(warn, "HMPID DigitToClusterSpec::initFileIn() : Did not find {} tree in file ", secondTree.c_str());
    return;
    std::exit(0);
  }

}


std::vector<o2::dataformats::MatchInfoHMP>* readTrack(int eventID, int trackID, int& pdg) {
  TFile* fMatch = new TFile("o2match_hmp.root");
  TTree* tMatch = (TTree*)fKine->Get("matchHMP");
	if(!tMatch) tMatch = (TTree*)fKine->Get("o2hmp");
  std::vector<o2::dataformats::MatchInfoHMP>* matchArr = nullptr;
  tMatch->SetBranchAddress("HMPMatchInfo", &matchArr);
  tMatch->GetEntry(eventID);
	
	Printf("matchArr size %d"  , matchArr->size());
	int cnt = 0;
  for(const auto& c : *matchArr) {	

		float x, y; int q, nph;
		c.getHMPIDmip(x, y, q, nph); 
	  Printf("i %d MIPindex %d" , cnt++, q);
  }

	return matchArr;
}

TTree* initializeTree(std::vector<o2::dataformats::MatchInfoHMP>*& matchArr) {
    TFile* fMatch = TFile::Open("o2match_hmp.root");
    if (!fMatch || fMatch->IsZombie()) {
        Printf("Error opening file");
        return nullptr;
    }
  
    TTree* tMatch = (TTree*)fMatch->Get("matchHMP");
    if(!tMatch) tMatch = (TTree*)fMatch->Get("o2hmp");
    if(!tMatch) {
        Printf("Error accessing TTree");
        fMatch->Close();
        delete fMatch;
        return nullptr;
    }

    tMatch->SetBranchAddress("HMPMatchInfo", &matchArr);
    tMatch->GetEntry(0);

    return tMatch;
}


// eventId = eventID to be searched for;
// startIndex : index of where matchArr is to be searched
// newStartIndex startIndex for next event
std::vector<o2::dataformats::MatchInfoHMP>* readTrackWithCut(TTree* tMatch, std::vector<o2::dataformats::MatchInfoHMP>* matchArr, int eventID, int& startIndex) {
    if(!tMatch) {
        Printf("TTree not initialized");
        return nullptr;
    }

    // Prepare to store filtered matches
    std::vector<o2::dataformats::MatchInfoHMP>* filteredMatches = new std::vector<o2::dataformats::MatchInfoHMP>;



    // tracks should be stored in "time" --> when we find our event we can then switch this condition "of" when the event changes:
    bool found = false;


    if((*matchArr)[startIndex].getEventNumber() != eventID) {Printf("This shouldnt happen");}
    else found = true;

    for(int i = startIndex; i < matchArr.size(); i++) {
        const auto& track = (*matchArr)[i];
        if(track.getEventNumber() != eventID) {
          startIndex = i; // new startIndex for next event
          break;
        } else { 
          filteredMatches->push_back(track);
        }
    }

    return filteredMatches;
}


TFile* fClu = new TFile("hmpidclusters.root");
TTree* tClu = (TTree*)fClu->Get("o2sim");
if(!tClu) tClu = (TTree*)fClu->Get("o2hmp");
std::vector<Cluster>* cluArr = nullptr;
tClu->SetBranchAddress("HMPIDclusters",&cluArr);

std::vector<Cluster>* readClu(int eventID, int trackID, int& pdg) {

  tClu->GetEntry(eventID);
	Printf("mcArr size %d"  , cluArr->size());

	int cnt = 0;
  for(const auto& c : *cluArr) {	
	  Printf("i %d ye %f" , cnt++, c.ye());
  }

	return cluArr;
}


std::vector<o2::MCTrack>* readTOF(int eventID, int trackID, int& pdg) {
  TFile* fKine = new TFile("o2sim_Kine.root");
  TTree* tKine = (TTree*)fKine->Get("o2sim");
  std::vector<o2::MCTrack>* mcArr = nullptr;
  tKine->SetBranchAddress("MCTrack", &mcArr);
  tKine->GetEntry(eventID);

	const auto& track = mcArr->at(trackID);



.

  for (int i = 0; i < mcArr->size(); ++i) {
    const auto& mcTrack = (*mcArr)[i];
    if (i == trackID) {
      Printf("Particle %d: pdg = %d, pT = %f, px = %f, py = %f, pz = %f, vx = %f, vy = %f, vz = %f", i, mcTrack.GetPdgCode(), TMath::Abs(mcTrack.GetStartVertexMomentumX() * mcTrack.GetStartVertexMomentumX() + mcTrack.GetStartVertexMomentumY() * mcTrack.GetStartVertexMomentumY()), mcTrack.GetStartVertexMomentumX(), mcTrack.GetStartVertexMomentumY(), mcTrack.GetStartVertexMomentumZ(), mcTrack.GetStartVertexCoordinatesX(), mcTrack.GetStartVertexCoordinatesY(), mcTrack.GetStartVertexCoordinatesZ());
    }
  }
	Printf("Particle  pdg = %d"  ,track.GetPdgCode()); 

	return mcArr;
}

void readTreeEntries() {
    // Open the ROOT file


    /*
    auto matchFile = std::make_unique<TFile>("o2match_hmp.root");
    auto clusterFile = std::make_unique<TFile>("hmpidclusters.root");
    auto mcFile = std::make_unique<TFile>("o2sim_Kine.root");

    std::unique_ptr<TTree> matchTree, clusTree, kineTRee; ///< input tree

    std::vector<Clusters>* mClustersFromFile;
    std::vector<o2::dataformats::MatchInfoHMP>* mTracksFromFile;
    std::vector<o2::MCTrack>* mCarloFromFile;
    initFileIn(matchFile, matchTree,  "matchHMP",  "matchHMP", "HMPMatchInfo", "HMPMatchInfo");
    /*if (mUseMC) {
        mTree->SetBranchAddress("MatchHMPMCTruth", &mLabelHMPPtr);
    } 

    initFileIn(clusterFile, clusTree,  "o2hmp",  "o2sim", "HMPIDClusters", "HMPIDclusters");*/



    /*void initFileIn(std::unique_ptr<TFile>& mFile, std::unique_ptr<TTree>& mTree, const std::string& firstTree, const std::string& secondTree, const std::string& firstBranch, const std::string& secondBranch, std::vector<T>*& dataFromFile)

    * / 



  if (clusTree->GetBranchStatus("HMPIDClusters") == 1) {
    clusTree->SetBranchAddress("HMPIDClusters", &mClustersFromFile);
  } else if (clusTree->GetBranchStatus("HMPIDclusters") == 1) {
    clusTree->SetBranchAddress("HMPIDclusters", &mClustersFromFile);
  } else {
  Printf("HMPID DigitToClusterSpec::initFileIn() : Error in branches!");
    return;
    std::exit(0);
  }


  if (matchTree->GetBranchStatus("HMPMatchInfo") == 1) {
    matchTree->SetBranchAddress("HMPMatchInfo", &mTracksFromFile);
  } else if (matchTree->GetBranchStatus("HMPMatchInfo") == 1) {
    matchTree->SetBranchAddress("HMPMatchInfo", &mTracksFromFile);
  } else {
    Printf("HMPID DigitToClusterSpec::initFileIn() : Error in branches!");
    return;
    std::exit(0);
  } */
 


	int i;
  auto trackVector = readTrack(0,0,i);
  auto clusterVector = readClu(0,0,i);
  std::vector<o2::MCTrack>* mcVector = readTOF(0,0,i);



  // 
  for() {
    
  }

	Printf("numTracks %d numClusters %d numMC %d", trackVector->size(), clusterVector->size(), mcVector->size());
	
  /*
  for(const auto& clu : *mClustersFromFile)
  {
	  const auto& dig = clu.dig(0);
		if(dig != nullptr) {

			const auto& trackId = dig->mTrackId;
			const auto& particlePdg = dig->mTrackId;

			int pdg;
			readTOF(0, trackId, pdg);
			Printf("trackId%d particlePdg %d pdg %d", trackId, particlePdg, pdg);
    }
	}*/ 
	


  //mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);

/* 
  matchTree->Print("toponly");
  clusTree->Print("toponly");
  //kineTRee->Print("toponly");
    // Loop over tree entries
    Long64_t nMatchEvents = matchTree->GetEntries();
    Long64_t nClusterEvents = clusTree->GetEntries();
    Long64_t nKineEvents = 5;//kineTRee->GetEntries();
    LOGP(info, " nMatchEvents {}, nClusterEvents {},  nKineEvents {}",nMatchEvents,nClusterEvents,nKineEvents);

    // loop over nKineTree see if we have HMP tracks and clusters for the given event


    //matchTree->GetEntry(0);  // This fills the above variables with data
    //clusTree->GetEntry(0);  // This fills the above variables with data


    LOGP(info, " numClusters = {}", mClustersFromFile->size());

	
  if (clusTree->GetReadEntry() + 1 >= clusTree->GetEntries()) {
    //mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
    //firstTrigger, lastTrigger = 0;
  } else {
    auto entry = clusTree->GetReadEntry() + 1;
    assert(entry < clusTree->GetEntries());
    clusTree->GetEntry(entry);

  }

  if (matchTree->GetReadEntry() + 1 >= matchTree->GetEntries()) {
    //mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
    //firstTrigger, lastTrigger = 0;
  } else {
    auto entry = matchTree->GetReadEntry() + 1;
    assert(entry < matchTree->GetEntries());
    matchTree->GetEntry(entry);
        
}*/
 //   LOGP(info, " nMatchEvents {}, nClusterEvents {},  nKineEvents {}",nMatchEvents,nClusterEvents,nKineEvents);





}
