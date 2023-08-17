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


TTree* initializeMatchTree(std::vector<o2::dataformats::MatchInfoHMP>* matchArr, int eventID, int trackID, int& pdg) {
  TFile* fMatch = new TFile("o2match_hmp.root");
  TTree* tMatch = (TTree*)fKine->Get("matchHMP");
	if(!tMatch) tMatch = (TTree*)fKine->Get("o2hmp");
  std::vector<o2::dataformats::MatchInfoHMP>* matchArr = nullptr;
  tMatch->SetBranchAddress("HMPMatchInfo", &matchArr);
  tMatch->GetEntry(0);
	return tMatch;
}


// eventId = eventID to be searched for;
// startIndex : index of where matchArr is to be searched
// newStartIndex startIndex for next event
std::vector<o2::dataformats::MatchInfoHMP>* readMatch(TTree* tMatch, std::vector<o2::dataformats::MatchInfoHMP>* matchArr, int eventID, int& startIndex) {
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



TTree* initializeClusterTree(std::vector<Cluster>*& cluArr) {
    TFile* fClu = TFile::Open("hmpidclusters.root");
    if (!fClu || fClu->IsZombie()) {
        Printf("Error opening file");
        return nullptr;
    }

    TTree* tClu = (TTree*)fClu->Get("o2sim");
    if(!tClu) tClu = (TTree*)fClu->Get("o2hmp");
    if(!tClu) {
        Printf("Error accessing TTree");
        fClu->Close();
        delete fClu;
        return nullptr;
    }

    tClu->SetBranchAddress("HMPIDclusters", &cluArr);
    tClu->GetEntry(0);
    return tClu;
}

TTree* initializeMCTree(std::vector<o2::MCTrack>*& mcArr) {
    TFile* fKine = TFile::Open("o2sim_Kine.root");
    if (!fKine || fKine->IsZombie()) {
        Printf("Error opening file");
        return nullptr;
    }

    TTree* tKine = (TTree*)fKine->Get("o2sim");
    if(!tKine) {
        Printf("Error accessing TTree");
        fKine->Close();
        delete fTrack;
        return nullptr;
    }

    tKine->SetBranchAddress("MCTrack", &mcArr);
    tKine->GetEntry(0);
    return tKine;
}

std::vector<o2::MCTrack>* readMC(std::vector<o2::MCTrack>*& mcArr, TTree* tKine, int eventId) {
  tKine->GetEntry(eventID);
	return mcArr;
}


// for the given eventID; read trackID 
o2::MCTrack getMCEntry(std::vector<o2::MCTrack>*& mcArr, int trackId) {
	const auto& track = mcArr->at(trackID);.
  for (int i = 0; i < mcArr->size(); ++i) {
    const auto& mcTrack = (*mcArr)[i];
    if (i == trackID) {
      Printf("Particle %d: pdg = %d, pT = %f, px = %f, py = %f, pz = %f, vx = %f, vy = %f, vz = %f", i, mcTrack.GetPdgCode(), TMath::Abs(mcTrack.GetStartVertexMomentumX() * mcTrack.GetStartVertexMomentumX() + mcTrack.GetStartVertexMomentumY() * mcTrack.GetStartVertexMomentumY()), mcTrack.GetStartVertexMomentumX(), mcTrack.GetStartVertexMomentumY(), mcTrack.GetStartVertexMomentumZ(), mcTrack.GetStartVertexCoordinatesX(), mcTrack.GetStartVertexCoordinatesY(), mcTrack.GetStartVertexCoordinatesZ());
      return mcTrack;
    }
  } 
  return nullptr;
}

void readTreeEntries() {
    // Open the ROOT file

	int i;
  auto trackVector = readTrack(0,0,i);
  auto clusterVector = readClu(0,0,i);
  std::vector<o2::MCTrack>* mcVector = readMC(0,0,i);
  for() {
    
  }

	Printf("numTracks %d numClusters %d numMC %d", trackVector->size(), clusterVector->size(), mcVector->size());

}
