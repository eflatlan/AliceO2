#ifndef HMPID_DATA_READER_H
#define HMPID_DATA_READER_H


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


class HmpidDataReader {
private:
    const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
    std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
    const float mass_Pion_sq = mass_Pion*mass_Pion, mass_Kaon_sq = mass_Kaon*mass_Kaon,  mass_Proton_sq = mass_Proton*mass_Proton;

public:
    HmpidDataReader() {
        // Constructor to initialize any member variables or resources
    }

    ~HmpidDataReader() {
        // Destructor to clean up any allocated resources
    }

    TTree* initializeMatchTree(std::vector<o2::dataformats::MatchInfoHMP>* matchArr, int eventID, int trackID, int& pdg);
    std::vector<o2::dataformats::MatchInfoHMP>* readMatch(TTree* tMatch, std::vector<o2::dataformats::MatchInfoHMP>* matchArr, int eventID, int& startIndex);
    TTree* initializeClusterTree(std::vector<Cluster>*& cluArr, std::vector<Trigger>*& trigArr);
    TTree* initializeMCTree(std::vector<o2::MCTrack>*& mcArr);
    std::vector<o2::MCTrack>* readMC(std::vector<o2::MCTrack>*& mcArr, TTree* tKine, int eventId);
    const o2::MCTrack* getMCEntry(std::vector<o2::MCTrack>* mcArr, int trackID);
    void readTreeEntries();
    std::array<float, 3> calcCherenkovHyp(float p, float n);
};

TTree* HmpidDataReader::initializeMatchTree(std::vector<o2::dataformats::MatchInfoHMP>* matchArr, int eventID, int trackID, int& pdg) {
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
std::vector<o2::dataformats::MatchInfoHMP>* HmpidDataReader::readMatch(TTree* tMatch, std::vector<o2::dataformats::MatchInfoHMP>* matchArr, int eventID, int& startIndex) {
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


    Printf("readMatch : event %d startIndex %d", eventID, startIndex);

    for(int i = startIndex; i < matchArr.size(); i++) {
        const auto& track = (*matchArr)[i];
        if(track.getEventNumber() != eventID) {
          startIndex = i; // new startIndex for next event
          Printf("readMatch : eventID changed - %d; end of loop ", track.getEventNumber());

          break;
        } else { 
          filteredMatches->push_back(track);
        }
    }

    Printf("readMatch : new startIndex %d", startIndex);

    return filteredMatches;
}



TTree* HmpidDataReader::initializeClusterTree(std::vector<Cluster>*& cluArr, std::vector<Trigger>*& trigArr) {
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
    tClu->SetBranchAddress("HMPIDInteractionRecords", &trigArr);

    tClu->GetEntry(0);
    return tClu;
}

TTree* HmpidDataReader::initializeMCTree(std::vector<o2::MCTrack>*& mcArr) {
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

std::vector<o2::MCTrack>* HmpidDataReader::readMC(std::vector<o2::MCTrack>*& mcArr, TTree* tKine, int eventId) {
  if (!fKine || fKine->IsZombie()) {
      Printf("Error opening file");
      return nullptr;
  }

  if(eventId < tKine->GetEentries()) {
    tKine->GetEntry(eventID);
    Printf("readMC at entry %d", eventID);
    return mcArr;
  } else {
    Printf("eventId > tKine->GetEntries()");
    return nullptr;
  }
	
}


// for the given eventID; read trackID 
const o2::MCTrack* HmpidDataReader::getMCEntry(std::vector<o2::MCTrack>* mcArr, int trackID) {

  if(trackID < 0 || trackID >= mcArr->size()){
    return nullptr;
  }
  else {
    //const auto& track = mcArr->at(trackID);
    for (int i = 0; i < mcArr->size(); ++i) {
      const auto& mcTrack = (*mcArr)[i];
      if (i == trackID) {
        Printf("Particle %d: pdg = %d, pT = %f, px = %f, py = %f, pz = %f, vx = %f, vy = %f, vz = %f", i, mcTrack.GetPdgCode(), TMath::Abs(mcTrack.GetStartVertexMomentumX() * mcTrack.GetStartVertexMomentumX() + mcTrack.GetStartVertexMomentumY() * mcTrack.GetStartVertexMomentumY()), mcTrack.GetStartVertexMomentumX(), mcTrack.GetStartVertexMomentumY(), mcTrack.GetStartVertexMomentumZ(), mcTrack.GetStartVertexCoordinatesX(), mcTrack.GetStartVertexCoordinatesY(), mcTrack.GetStartVertexCoordinatesZ());
        return &mcTrack;
      }
    }
  }
  return nullptr;
}

void HmpidDataReader::readTreeEntries() {
    // Open the ROOT file

	int i;
  auto trackVector = readTrack(0,0,i);
  auto clusterVector = readClu(0,0,i);
  std::vector<o2::MCTrack>* mcVector = readMC(0,0,i);
  for() {
    
  }

	Printf("numTracks %d numClusters %d numMC %d", trackVector->size(), clusterVector->size(), mcVector->size());

}


std::array<float, 3> HmpidDataReader::calcCherenkovHyp(float p, float n)
{
  const float p_sq = p*p;
  const float cos_ckov_denom = p*refIndexFreon;
  const auto cos_ckov_Pion = static_cast<float>(TMath::Sqrt(p_sq + mass_Pion_sq)/(cos_ckov_denom)); // n = refIndexFreon 1.289 later make it random?

  const float cos_ckov_Kaon = static_cast<float>(TMath::Sqrt(p_sq + mass_Kaon_sq)/(cos_ckov_denom)); 
  const float cos_ckov_Proton = static_cast<float>(TMath::Sqrt(p_sq + mass_Proton_sq)/(cos_ckov_denom));
  

  const float ckovAnglePion = static_cast<float>( TMath::ACos(cos_ckov_Pion)); 
  const float ckovAngleKaon = static_cast<float>(TMath::ACos(cos_ckov_Kaon)); 
  const float ckovAngleProton = static_cast<float>(TMath::ACos(cos_ckov_Proton)); 

  Printf("Pion %.3f Kaon %.3f Proton %.3f", ckovAnglePion, ckovAngleKaon, ckovAngleProton);

  return {ckovAnglePion, ckovAngleKaon, ckovAngleProton};
}

#endif // HMPID_DATA_READER_H
