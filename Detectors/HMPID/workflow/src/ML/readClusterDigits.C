#if !defined(__CLING__) || defined(__ROOTCLING__)
// Figures 4.5-4.6 
// ROOT header-files
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
#include <TROOT.h> // gRoot
#include <TGraph.h>

// O2 header-files
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "DataFormatsHMP/Cluster.h"

// C++ header files and libraries
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <ctime>    
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#endif

void SaveFolder(string inpFile);
void changeFont();

void readDigits(char* filename, int nEvent)
{
  TFile *fileDigits = TFile::Open(filename);
  // Cast to TTree*, get tree by key "o2sim"
  TTree *treeDigits = (TTree*)fileDigits->Get("o2sim"); 
  
  // Initialize an array of TH1F-pointers 
  // (TH1=1D histogram, F specifies Float-values )
  TH1F *hCharge[7];   	//  charge of digits
  TH1F *h_xCoord[7];	  //  number of digits per x-coordinate 
  TH1F *h_yCoord[7];    //  number of digits per y-coordinate

  // Initialize an array of TH2F-pointers 
  TH2F *hMap[7];	// 2d map of digits
      
  // Label the histograms
  for(int i=0; i<7; i++) {
    // define element number i in the pointer-array

    hMap[i] = new TH2F(Form("Digits Map %i",i),\
    Form("Digits Map %i",i), 160, 0, 159, 144, 0, 143);  
    hMap[i]->SetXTitle("pad X [cm]");
    hMap[i]->SetYTitle("pad Y [cm]"); 
   
    hCharge[i] = new TH1F(Form("Digits Charge %i",i),\
    Form("Digits Charge %i",i),2000, 100., 2100.);
    hCharge[i]->SetXTitle("Charge (ADC channel)");               
    hCharge[i]->SetYTitle("Entries");         

    h_xCoord[i] = new TH1F(Form("Digits X-location Histogram %i",i),\
    Form("Digits X-location Histogram %i",i),2000, 10., 159.);
    h_xCoord[i]->SetXTitle("X [cm]");               
    h_xCoord[i]->SetYTitle("Entries");  


    h_yCoord[i] = new TH1F(Form("Digits Y-location Histogram %i",i),\
    Form("Digits Y-location Histogram %i",i),2000, 10., 144.);
    h_yCoord[i]->SetXTitle("Y [cm]");               
    h_yCoord[i]->SetYTitle("Entries");    
     
  }
  
  changeFont();		         // apply custom canvas figure options
  //SaveFolder("digitChambers"); // specify folder to save files in

  // Define canvases for plotting the figures  
  TCanvas *c1 = new TCanvas("c1","c1",2000,1200); 	
  TCanvas *c2 = new TCanvas("c2","c2",2000,1200); 	
  TCanvas *c3 = new TCanvas("c3","c3",2000,1200); 	
  TCanvas *c4 = new TCanvas("c4","c4",2000,1200); 	

  // Define positions for the plots for the chambers in the canvases 
  Int_t pos[]= {9,8,6,5,4,2,1};
  
  o2::hmpid::Trigger *pTgr;  // pointer to Trigger-object
  o2::hmpid::Digit *pDig;    // pointer to Digit-object
  o2::hmpid::Digit *pDigEvt; // pointer to Digit-object
  o2::hmpid::Digit digit;    // declaration of digit-object
  
  std::vector<o2::hmpid::Digit> *digits = nullptr;    //  vector of digit-pointers
  std::vector<o2::hmpid::Digit> oneEventDigits;       //  vector of Digit-objects
  std::vector<o2::hmpid::Trigger> *trigger = nullptr; //  vector of trigger-pointers
  
  // specify branches of tree to read from 
  // (const char* branchname, void* adress)
  treeDigits->SetBranchAddress("HMPDigit",&digits); 
  treeDigits->SetBranchAddress("InteractionRecords",&trigger);

  // Number of entries in tree
  Printf("tree entries = %i", treeDigits->GetEntries()); 
  
  treeDigits->GetEntry(0);

  Printf("digit size = %i", digits->size());
  
  // initialize padhit x and y coordinates, 
  // and module (i.e.) chamber-number
  int padChX = 0, padChY = 0, module = 0;

  // Loop through digits in file
  for(unsigned int j = 0; j < digits->size(); j++) {
      
    pDig = (o2::hmpid::Digit*)&digits->at(j); 
    // digits->at(j) = use arrow-operator since it is an array of pointers
    //&digits = create pointer pDig to adress of vector-element j in digits 
    // (o2::hmpid::Digit*) = cast to Digit*
    
    // Call member-function pad2Absolute in Digit-class
    // static void pad2Absolute(uint32_t pad, int* Module, int* x, int* y);
    // No need to be called from an object since the member function is static
    o2::hmpid::Digit::pad2Absolute(pDig->getPadID(), &module, &padChX, &padChY);    
    
    // Fill histograms, with the module specifying which chamber, 
    //(= pad in the canvas), should be filled 
    hCharge[module]->Fill(pDig->getQ()); // getQ gets charge
    h_xCoord[module]->Fill(padChX);	     // fill x-coordinate 
    h_yCoord[module]->Fill(padChY);	     // fill y-coordinate
  }
       
  Printf("trigger size from digits = %i", trigger->size()); 
   
	
  // Loop through triggers 
  for(int i = 0; i < trigger->size(); i++) {
    
    oneEventDigits.clear();	// empty vector 
    pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
    
    for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      // get all digits corresponding to one trigger, and assemble them 
      // in a vector
      digit = (o2::hmpid::Digit)digits->at(j);
      oneEventDigits.push_back(digit);
            
    }
   
    // get 2D-map for the event-number specified as input   
    if(i==nEvent) {
       
      for(unsigned int k = 0; k < oneEventDigits.size(); k++) {
      
        pDigEvt = (o2::hmpid::Digit*)&oneEventDigits.at(k);
                     
        o2::hmpid::Digit::pad2Absolute(pDigEvt->getPadID(), &module, &padChX, &padChY);
         
        hMap[module]->Fill(padChX, padChY, pDigEvt->getQ()); 
      LOGP(info, "dig In event ; {}", k);			                       
      }              
    }// end if 
  }
  // Divide canvases into 3x3     
  c1->Divide(3,3);
  c2->Divide(3,3); 
  c3->Divide(3,3);  	
  c4->Divide(3,3);  
	


 
 for(int iCh = 0; iCh<7; iCh++){
  
    c2->cd(pos[iCh]); 
    hCharge[iCh]->SetStats(false);
    hCharge[iCh]->Draw();
    
    c1->cd(pos[iCh]); 
    hMap[iCh]->SetStats(false);
    hMap[iCh]->Draw("Colz");
     
    c3->cd(pos[iCh]); 
    h_xCoord[iCh]->SetStats(false);
    h_xCoord[iCh]->Draw();
   
    c4->cd(pos[iCh]); 
    h_yCoord[iCh]->SetStats(false);
    h_yCoord[iCh]->Draw();

  }      
  c1->SaveAs("digitmap.eps");
  c2->SaveAs("digitcharge.eps");
  c3->SaveAs("digitX.eps");
  c4->SaveAs("digitY.eps");      
}
//************************************************************************************************************************************************************************************************************************************************************************************
void readClusters(char* filename = "hmpidclusters.root", int nEvent = 1)
{ 
  TFile *fileClusters = TFile::Open(filename);
  
  // get tree from root file with key "o2hmp"
  TTree *treeClusters = (TTree*)fileClusters->Get("o2hmp");
  
  // Initialize an array of TH1F-pointers 
  //(TH1=1D histogram, F specifies Float-values )
  TH1F *hCharge[7]; 	     // charge of clusters
  TH1F *hSize[7]; 	     // size of clusters
  TH1F *h_xCoord[7];	     // number of custers per x-coordinate 
  TH1F *h_yCoord[7];         //  number of clusters per y-coordinate

  TH2F *hMap[7]; 	     // 2d histogram/ map of clusters

  // Define positions for the plots for the chambers in the canvases 
  Int_t pos[]= {9,8,6,5,4,2,1};
   

  // Label the histograms
  for(int i=0; i<7; i++) {

    hMap[i] = new TH2F(Form("Clusters Map %i",i),\
    Form("Clusters Map %i",i), 160, 0, 159, 144, 0, 143);
    hMap[i]->SetXTitle("X (cm)");
    hMap[i]->SetYTitle("Y (cm) "); 
     
    hCharge[i] = new TH1F(Form("Clusters Charge %i",i),\
    Form("Clusters Charge %i",i),2000, 100., 2100.);
    hCharge[i]->SetXTitle("Charge (ADC channel)");               
    hCharge[i]->SetYTitle("Entries ");
     
    hSize[i] = new TH1F(Form("Cluster Size chmaber%i",i),\
    Form("Cluster Size %i",i),20, 0., 20.);
    hSize[i]->SetXTitle("Cluster size");               
    hSize[i]->SetYTitle("Entries ");               
     
    h_xCoord[i] = new TH1F(Form("Cluster X-location Histogram %i",i),\
    Form("Cluster X-location Histogram %i",i),2000, 10., 159.);
    h_xCoord[i]->SetXTitle("X [cm]");               
    h_xCoord[i]->SetYTitle("Entries ");  

    h_yCoord[i] = new TH1F(Form("Cluster Y-location Histogram %i",i),\
    Form("Cluster Y-location Histogram %i",i),2000, 10., 144.);
    h_yCoord[i]->SetXTitle("Y [cm]");               
    h_yCoord[i]->SetYTitle("Entries ");                 
  }
     
  changeFont();			     // specify folder to save files in
  //SaveFolder("clusterChambers");   // apply custom canvas figure options

  o2::hmpid::Trigger *pTgr;       // pointer to Trigger-object
  o2::hmpid::Cluster *pClu;       // pointer to cluster-object
  o2::hmpid::Cluster *pCluEvt;    // pointer to cluster-object
  o2::hmpid::Cluster cluster;     // declaration to cluster-object
  
  std::vector<o2::hmpid::Cluster> *clusters = nullptr; //  vector of cluster-pointers
  std::vector<o2::hmpid::Cluster> oneEventClusters;    //  vector of cluster-objects
  std::vector<o2::hmpid::Trigger> *trigger = nullptr;  //  vector of Trigger-pointers
  
  // specify branches of tree to read from 
  treeClusters->SetBranchAddress("HMPIDClusters",&clusters);

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
  
  Printf("tree entries = %i", treeClusters->GetEntries());
  
  treeClusters->GetEntry(0);

  Printf("clusters size = %i", clusters->size());

  // initialize module (i.e. chamber-number [0..6])
  int module = 0;
  
  // loop through clusters in file
  for(unsigned int j = 0; j < clusters->size(); j++) {
      
    pClu = (o2::hmpid::Cluster*)&clusters->at(j);
	
    module = pClu->ch(); // get module (chamber-number [0..6]) for the given cluster
                             
    // fill 1d histograms
    hCharge[module]->Fill(pClu->q());  // get charge for the cluster     
    hSize[module]->Fill(pClu->size()); // get size of cluster 
    h_xCoord[module]->Fill(pClu->x()); // cluster x position in LRS
    h_yCoord[module]->Fill(pClu->y()); // cluster y position in LRS

  }
        
  Printf("trigger size from clusters = %i", trigger->size()); 
   
  // loop through triggers
  for(int i = 0; i < trigger->size(); i++) {

    oneEventClusters.clear();	// empty vector

    pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
    
    for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      // get all clusters corresponding to one trigger, and assemble them 
      // in a vector
      cluster = (o2::hmpid::Cluster)clusters->at(j);
      oneEventClusters.push_back(cluster);
      
     }
    
    // get 2D-map for the event-number specified as input     
    if(i==nEvent) {
       
      for(unsigned int k = 0; k < oneEventClusters.size(); k++) {
      
      pCluEvt = (o2::hmpid::Cluster*)&oneEventClusters.at(k);
                     
      module = pCluEvt->ch();
       
      hMap[module]->Fill(pCluEvt->x(), pCluEvt->y());
      LOGP(info, "Cluster In event ; {}", k);			     
      }              
    }//end if
  } 

  // define canvas-names 
  TCanvas *c1 = new TCanvas("c1","c1",2000,1200); 	
  TCanvas *c2 = new TCanvas("c2","c2",2000,1200); 	
  TCanvas *c3 = new TCanvas("c3","c3",2000,1200); 	
  TCanvas *c4 = new TCanvas("c4","c4",2000,1200); 	
  TCanvas *c5 = new TCanvas("c5","c5",2000,1200); 	
   
  // divide canvases into 3x3 pads
  c1->Divide(3,3);
  c2->Divide(3,3);
  c3->Divide(3,3);
  c4->Divide(3,3);
  c5->Divide(3,3);

	
  // draw canvases and save files
  for(int iCh = 0; iCh<7; iCh++){
    c1->cd(pos[iCh]); 
    hMap[iCh]->SetStats(false);
    hMap[iCh]->SetMarkerStyle(3);
    hMap[iCh]->Draw();
     
    c2->cd(pos[iCh]); 
    hMap[iCh]->SetStats(false);
    hCharge[iCh]->Draw();
     
    c3->cd(pos[iCh]); 
    hSize[iCh]->SetStats(false);
    hSize[iCh]->Draw();

    c4->cd(pos[iCh]); 
    h_xCoord[iCh]->SetStats(false);
    h_xCoord[iCh]->Draw();

    c5->cd(pos[iCh]); 
    h_yCoord[iCh]->SetStats(false);
    h_yCoord[iCh]->Draw();
    
  }       
  c1->SaveAs("clustermap.eps");
  c2->SaveAs("clustercharge.eps");
  c3->SaveAs("clustersize.eps");
  c4->SaveAs("clusterX.eps");
  c5->SaveAs("clusterY.eps");
}    

//**********************************************************************************************


// save files in specified folder    
void SaveFolder(string inpFile)
{
	 char* createdFolder;
	 auto time = std::chrono::system_clock::now();
	 std::time_t time_t = std::chrono::system_clock::to_time_t(time);
	 auto c_time = std::ctime(&time_t);
	
	 // Allocate memory for char* createdFolder 
	 int  pathLen = strlen(c_time);
	 pathLen = strlen(gSystem->pwd())+strlen(inpFile.c_str());
	 int numOfSigns = 2; // allocate for - and /
	 createdFolder = (char *)calloc(pathLen+numOfSigns+1, sizeof(char));

	 // Copy Base-directory into empty createdFolder 
	 strcpy(createdFolder,gSystem->pwd());
	 // add hyphen to the Folder-name  
	 strcat(createdFolder,"/"); 
	 // Add name of read file to Folder-name
	 strcat(createdFolder,inpFile.c_str()); 
	 // add - to the Folder-name
	 strcat(createdFolder,"-");
	 // append current time to the Folder-name
	 strcat(createdFolder,c_time); 

	 // Make new Directory from the newly 
	 //created char* createdFolder
	 gSystem->MakeDirectory(createdFolder);  
		
	 // Move to directory such that new
	 // canvas-files will be saved here  
	 gSystem->ChangeDirectory(createdFolder);  
		
}


// apply custom-made canvas options 
void changeFont()
{
  TStyle* canvasStyle = new TStyle("canvasStyle","Canvas Root Styles");
  canvasStyle->SetPalette(1,0);
  canvasStyle->SetTitleSize(0.085,"xy");    // size of axis title font
  canvasStyle->SetTitleFont(22,"xz");       // font option
  canvasStyle->SetTitleFontSize(0.1);       // size of canvas-title
  canvasStyle->SetTitleOffset(.825,"y");    //  y-axis title-offset from axis
  canvasStyle->SetTitleOffset(1,"z"); 
  canvasStyle->SetTitleOffset(.95,"x");
  canvasStyle->SetTitleX(.25);              // set canvas-title position
  // labels 
  canvasStyle->SetLabelOffset(0.005,"y"); 
  canvasStyle->SetLabelFont(22,"xyz");
  canvasStyle->SetLabelSize(0.085,"xyz");   // size of axis value font
  // canvas 
  canvasStyle->SetCanvasColor(0); 
  canvasStyle->SetCanvasBorderMode(0);
  canvasStyle->SetCanvasBorderSize(0);
  //set margins
  canvasStyle->SetPadBottomMargin(0.18);
  canvasStyle->SetPadTopMargin(0.05);
  canvasStyle->SetPadLeftMargin(0.13);
  canvasStyle->SetPadRightMargin(0.02);
  gROOT->SetStyle("canvasStyle"); 
}

