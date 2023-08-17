#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TList.h>
#include <TLine.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TGraph.h>
#include <fstream>
#include <vector>
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "DataFormatsHMP/Cluster.h"
#endif


using Cluster = o2::hmpid::Cluster;


void readDigits(char* filename, int nEvent)
{
  TFile *fileDigits = TFile::Open(filename);
    std::map<int, int> pdgCounts;  // Map to store PDG values and their counts

  TTree *treeDigits = (TTree*)fileDigits->Get("o2sim");
  
  TH1F *hCharge[7];
  TH2F *hMap[7];
      
    TH2F *hParticlePdgMapAbove[7];
    TH2F *hParticlePdgMapBelow[7];


TH2F* hParticlePdgElectron[7], *hParticlePdgPion[7], *hParticlePdgKaon[7],  *hParticlePdgProton[7], *hParticlePdgPhoton[7];


//11 111 321 2212 22/50000050
for (int i = 0; i < 7; i++) {
    // Electron histograms
    hParticlePdgElectron[i] = new TH2F(Form("Electron_ParticlePdg_Chamber%i", i),
                                       Form("Electron PDG Map Chamber%i", i),
                                       160, 0, 159, 144, 0, 143);
    hParticlePdgElectron[i]->SetXTitle("pad X");
    hParticlePdgElectron[i]->SetYTitle("pad Y");
    hParticlePdgElectron[i]->SetMarkerColor(kBlue); // Example color

    // Pion histograms
    hParticlePdgPion[i] = new TH2F(Form("Pion_ParticlePdg_Chamber%i", i),
                                   Form("Pion PDG Map Chamber%i", i),
                                   160, 0, 159, 144, 0, 143);
    hParticlePdgPion[i]->SetXTitle("pad X");
    hParticlePdgPion[i]->SetYTitle("pad Y");
    hParticlePdgPion[i]->SetMarkerColor(kRed); // Example color
    //hParticlePdgPion[i]->SetMarkerStyle(3); // Example color
    // Kaon histograms
    hParticlePdgKaon[i] = new TH2F(Form("Kaon_ParticlePdg_Chamber%i", i),
                                   Form("Kaon PDG Map Chamber%i", i),
                                   160, 0, 159, 144, 0, 143);
    hParticlePdgKaon[i]->SetXTitle("pad X");
    hParticlePdgKaon[i]->SetYTitle("pad Y");
    hParticlePdgKaon[i]->SetMarkerColor(kGreen); // Example color
    //hParticlePdgKaon[i]->SetMarkerStyle(3); // Example color
    // Proton histograms
    hParticlePdgProton[i] = new TH2F(Form("Proton_ParticlePdg_Chamber%i", i),
                                     Form("Proton PDG Map Chamber%i", i),
                                     160, 0, 159, 144, 0, 143);
    hParticlePdgProton[i]->SetXTitle("pad X");
    hParticlePdgProton[i]->SetYTitle("pad Y");
    hParticlePdgProton[i]->SetMarkerColor(kOrange); // Example color
    hParticlePdgProton[i]->SetMarkerStyle(3); // Example color
    // Photon histograms
    hParticlePdgPhoton[i] = new TH2F(Form("Photon_ParticlePdg_Chamber%i", i),
                                     Form("Photon PDG Map Chamber%i", i),
                                     160, 0, 159, 144, 0, 143);
    hParticlePdgPhoton[i]->SetXTitle("pad X");
    hParticlePdgPhoton[i]->SetYTitle("pad Y");
    hParticlePdgPhoton[i]->SetMarkerColor(kBlack); // Example color
}
// Add a new histogram to collect PDG values
TH1F* hPDG = new TH1F("hPDG", "Particle PDG values", 5000, -2500, 2500);
hPDG->SetXTitle("PDG Value");
hPDG->SetYTitle("Entries");


    for (int i = 0; i < 7; i++) {
        hParticlePdgMapAbove[i] = new TH2F(Form("Digits ParticlePdg Map Above 10000 chamber%i", i), 
                                          Form("Digits ParticlePdg Map Above 10000 chamber%i", i), 
                                          160, 0, 159, 144, 0, 143);
        hParticlePdgMapAbove[i]->SetXTitle("pad X");
        hParticlePdgMapAbove[i]->SetYTitle("pad Y");
        hParticlePdgMapAbove[i]->SetMarkerColor(kRed);

        hParticlePdgMapBelow[i] = new TH2F(Form("Digits ParticlePdg Map Below 10000 chamber%i", i), 
                                          Form("Digits ParticlePdg Map Below 10000 chamber%i", i), 
                                          160, 0, 159, 144, 0, 143);
        hParticlePdgMapBelow[i]->SetXTitle("pad X");
        hParticlePdgMapBelow[i]->SetYTitle("pad Y");
        hParticlePdgMapBelow[i]->SetMarkerColor(kGreen);
    }

  for(int i=0; i<7; i++) {
     hMap[i] = new TH2F(Form("Digits Map chmaber%i",i),Form("Digits Map chmaber%i",i), 160, 0, 159, 144, 0, 143);
     hMap[i]->SetXTitle("pad X");
     hMap[i]->SetYTitle("pad Y"); 
     
     hCharge[i] = new TH1F(Form("Digits Charge chamber%i",i),Form("Digits Charge chamber%i",i),2000, 100., 2100.);
     hCharge[i]->SetXTitle("Charge (ADC channel)");               
     hCharge[i]->SetYTitle("Entries");               
   }
   
  TCanvas *c1 = new TCanvas("c1","c1",1000,800); 
  TCanvas *c2 = new TCanvas("c2","c2",1000,800); 
   
  Int_t pos[]= {9,8,6,5,4,2,1};
  
  o2::hmpid::Trigger *pTgr;
  o2::hmpid::Digit *pDig;
  o2::hmpid::Digit *pDigEvt;
  o2::hmpid::Digit digit;
  
  std::vector<o2::hmpid::Digit> *digits = nullptr;
  std::vector<o2::hmpid::Digit> oneEventDigits;
  std::vector<o2::hmpid::Trigger> *trigger = nullptr;
  
  treeDigits->SetBranchAddress("HMPDigit",&digits);
  treeDigits->SetBranchAddress("InteractionRecords",&trigger);
  
  Printf("tree entries = %i", treeDigits->GetEntries());
  
  treeDigits->GetEntry(0);

  Printf("digit size = %i", digits->size());
  
  int padChX = 0, padChY = 0, module = 0;
    



//11 111 321 2212 22/50000050
  std::vector<int> particlePDGs;
  for(unsigned int j = 0; j < digits->size(); j++) {


    pDig = (o2::hmpid::Digit*)&digits->at(j);
    hPDG->Fill(pDig->mParticlePdg);
  pdgCounts[pDig->mParticlePdg]++;
    o2::hmpid::Digit::pad2Absolute(pDig->getPadID(), &module, &padChX, &padChY);    
                         
    hCharge[module]->Fill(pDig->getQ());

        if (pDig->mParticlePdg > 10000) {
            hParticlePdgMapAbove[module]->Fill(padChX, padChY);
        } else {
            hParticlePdgMapBelow[module]->Fill(padChX, padChY);
        }


    switch (TMath::Abs(pDig->mParticlePdg)) {
      case 11 : 
	hParticlePdgElectron[module]->Fill(padChX, padChY); 
	break;
      case 211 : 
	hParticlePdgPion[module]->Fill(padChX, padChY);
	break;
      case 321 : 
	hParticlePdgKaon[module]->Fill(padChX, padChY);
	break;
      case 2212 : 
	hParticlePdgProton[module]->Fill(padChX, padChY);
	break;
      case 50000050 : 
	hParticlePdgPhoton[module]->Fill(padChX, padChY);
	break;
    }
 
    bool isInPdg = false;
    for(const auto& pdg: particlePDGs) { if(pdg == pDig->mParticlePdg) {isInPdg = true;} }
    if(!isInPdg) {particlePDGs.push_back(pDig->mParticlePdg);}    

   }

  Printf("PDG code counts:");
  for (const auto& pair : pdgCounts) {
    Printf("PDG: %d, Count: %d", pair.first, pair.second);
  }
for(const auto& pdg: particlePDGs) { 
  Printf("pdg = %d", pdg); 
 }
       
  Printf("trigger size from digits = %i", trigger->size()); 
   
  for(int i = 0; i < trigger->size(); i++) {
    
    oneEventDigits.clear();
    pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
    
    for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      
      digit = (o2::hmpid::Digit)digits->at(j);
      oneEventDigits.push_back(digit);
      
     }
      
    if(i==nEvent) {
       
       for(unsigned int k = 0; k < oneEventDigits.size(); k++) {
      
       pDigEvt = (o2::hmpid::Digit*)&oneEventDigits.at(k);
                     
       o2::hmpid::Digit::pad2Absolute(pDigEvt->getPadID(), &module, &padChX, &padChY);
         
       hMap[module]->Fill(padChX, padChY, pDigEvt->getQ());
                    
      }              
     }
   }
       
   c1->Divide(3,3);
   c2->Divide(3,3);
   TCanvas* c3 = new TCanvas("c3", "Particle Pdg Map", 800, 600);
   c3->Divide(3, 3); 

   TCanvas* c4 = new TCanvas("c4", "Particle Pdg Map Below", 800, 600);
   c4->Divide(3, 3);

   TCanvas* c5 = new TCanvas("c5", "Particle Pdg Map", 800, 600);
   c5->Divide(3, 3);



TCanvas* c6 = new TCanvas("c6", "PDG Histogram", 800, 600);
	c6->cd();
	hPDG->Draw();

   for(int iCh = 0; iCh<7; iCh++){
     c1->cd(pos[iCh]); 
     hMap[iCh]->Draw("Colz");
     
     c2->cd(pos[iCh]); 
     hCharge[iCh]->Draw();

     c3->cd(pos[iCh]); 
     hParticlePdgMapBelow[iCh]->Draw();
     //hParticlePdgMapAbove[iCh]->Draw("Colz");

     c4->cd(pos[iCh]); 
     hParticlePdgMapAbove[iCh]->Draw();
     //hParticlePdgMapAbove[iCh]->Draw("Colz");

     c5->cd(pos[iCh]); 
     hParticlePdgElectron[iCh]->Draw();
     hParticlePdgPion[iCh]->Draw("same");
     hParticlePdgKaon[iCh]->Draw("same");	
     hParticlePdgProton[iCh]->Draw("same");
     hParticlePdgPhoton[iCh]->Draw("same");
   }      
}
//************************************************************************************************************************************************************************************************************************************************************************************
void readClusters(char* filename, int nEvent)
{ 
  TH1F *hCharge[7], *hMipCharge[7], *hSize[7];
  TH2F *hMap[7];
      
  for(int i=0; i<7; i++) {
     hMap[i] = new TH2F(Form("Clusters Map chamber%i",i),Form("Cluster Map chamber%i",i), 160, 0, 159, 144, 0, 143);
     hMap[i]->SetXTitle("X (cm)");
     hMap[i]->SetYTitle("Y (cm)"); 
     
     hCharge[i] = new TH1F(Form("Clusters Charge chamber%i",i),Form("Cluster Charge chamber%i",i),2000, 100., 2100.);
     hCharge[i]->SetXTitle("Charge (ADC channel)");               
     hCharge[i]->SetYTitle("Entries");
     
     hMipCharge[i] = new TH1F(Form("Mip Clusters Charge chamber%i",i),Form("Mip Cluster Charge chamber%i",i),50, 200., 2200.);
     hMipCharge[i]->SetXTitle("Charge (ADC channel)");               
     hMipCharge[i]->SetYTitle("Entries/40 ADC");
     hMipCharge[i]->SetLineColor(kBlack);
     
     hSize[i] = new TH1F(Form("Cluster Size chmaber%i",i),Form("Cluster Size chamber%i",i),20, 0., 20.);
     hSize[i]->SetXTitle("Cluster size");               
     hSize[i]->SetYTitle("Entries");               
                    
   }
   
  TCanvas *c1 = new TCanvas("c1","c1",1000,800); 
  TCanvas *c2 = new TCanvas("c2","c2",1000,800); 
  TCanvas *c3 = new TCanvas("c3","c3",1000,800); 
  TCanvas *c4 = new TCanvas("c4","c4",1000,800);
   
  Int_t pos[] = {9,8,6,5,4,2,1};
  
  Int_t nTotTriggers = 0;
  
  for(int k = 0; k<1; k++) {
    
    Printf("s", filename);
      
    TFile *fileClusters = TFile::Open(Form("%s", filename));
    //TFile *fileClusters = TFile::Open(Form("%s%03i.root", filename, k+1));
    
    TTree *treeClusters = (TTree*)fileClusters->Get("o2hmp");
    


if (!treeClusters) {
treeClusters = ((TTree *)fileClusters->Get("o2sim"));
}
    o2::hmpid::Trigger *pTgr;
    o2::hmpid::Cluster *pClu;
    o2::hmpid::Cluster *pCluEvt;
    o2::hmpid::Cluster cluster;
  
    std::vector<o2::hmpid::Cluster> *clusters = nullptr;
    std::vector<o2::hmpid::Cluster> oneEventClusters;
    std::vector<o2::hmpid::Trigger> *trigger = nullptr;
  
    if(!treeClusters) {
      Printf("Empty clusters");
      continue;
     }
     
    else       
      {
	  if ((treeClusters->GetBranchStatus("HMPIDClusters")) == 1) {
	    treeClusters->SetBranchAddress("HMPIDClusters", &clusters);
	  } else if ((treeClusters->GetBranchStatus("HMPIDclusters")) == 1) {
	    treeClusters->SetBranchAddress("HMPIDclusters", &clusters);
	  } else {
	    LOG(warn) << "HMPID HMPIDclusters::init() : Error in branches!"
		      << endl;
	    return;
	    std::exit(0);
	  }

       treeClusters->SetBranchAddress("InteractionRecords",&trigger);
      }

    Printf("tree entries = %i", treeClusters->GetEntries());
  
    treeClusters->GetEntry(0);

    Printf("clusters size = %i", clusters->size());
  
    int module = 0;
    
    for(unsigned int j = 0; j < clusters->size(); j++) {
      
      pClu = (o2::hmpid::Cluster*)&clusters->at(j);

      module = pClu->ch();
                             
      hCharge[module]->Fill(pClu->q());
      std::vector<Cluster::Topology> topologies = pClu->getTopologyVector();
      

      std::cout << "evnum" << pClu->getEventNumber()<< std::endl; 
      /*for (const Cluster::Topology& topo : topologies) {
              std::cout << "diffX: " << topo.diffX << ", "
              << "diffY: " << topo.diffY << ", "
              << "q: " << topo.q << ", "
              << "pdg: " << topo.pdg << ", "
              << "tid: " << topo.tid << ", "
              << "mid: " << topo.mid 
              << std::endl;
      }*/

    
      if(pClu->size() >=3 && pClu->size()<=7) hMipCharge[module]->Fill(pClu->q());
    
      hSize[module]->Fill(pClu->size());
   }
       
   nTotTriggers+=trigger->size();
   
   Printf("trigger size from clusters = %i", trigger->size()); 
   
   for(int i = 0; i < trigger->size(); i++) {
     oneEventClusters.clear();
     pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
    
     for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      
      cluster = (o2::hmpid::Cluster)clusters->at(j);
      oneEventClusters.push_back(cluster);
      
     }
      
      if(i==nEvent) {
       
       for(unsigned int k = 0; k < oneEventClusters.size(); k++) {
      
         pCluEvt = (o2::hmpid::Cluster*)&oneEventClusters.at(k);
                     
         module = pCluEvt->ch();
       
         hMap[module]->Fill(pCluEvt->x(), pCluEvt->y());
                    
       } 
      }            
     }
   }
   
   c1->Divide(3,3);
   c2->Divide(3,3);
   c3->Divide(3,3);
   c4->Divide(3,3);

   for(int iCh = 0; iCh<7; iCh++){
     c1->cd(pos[iCh]); 
     hMap[iCh]->SetMarkerStyle(3);
     hMap[iCh]->Draw();
     
     c2->cd(pos[iCh]); 
     hCharge[iCh]->Draw();
     
     c3->cd(pos[iCh]); 
     hMipCharge[iCh]->Draw();
     
     c4->cd(pos[iCh]); 
     hSize[iCh]->Draw();
     
   }  
   
  Printf("Nuomber of triggers = %i",nTotTriggers);       
}    
