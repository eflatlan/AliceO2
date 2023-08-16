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
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"

using Digit = o2::hmpid::Digit;
using Cluster = o2::hmpid::Cluster;


void readDigits(char* filename, int nEvent) 
{
    TFile* f = TFile::Open(filename);
    if (!f) {
        std::cout << "Error: Cannot open file!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)f->Get("o2hmp");
    


		if (!tree) {
		tree = ((TTree *)f->Get("o2sim"));
		}
    std::vector<Digit>* digits = nullptr;
    tree->SetBranchAddress("HMPDigit", &digits);

    TList* list = new TList();

    TH2F *hParticlePdgMapAbove[7];
    TH2F *hParticlePdgMapBelow[7];
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

    int nEntries = tree->GetEntries();
    if (nEvent >= nEntries) {
        std::cout << "Error: Requested event is out of range!" << std::endl;
        return;
    }

    tree->GetEntry(nEvent);
    Printf("digit size = %i", digits->size());

    for (unsigned int k = 0; k < digits->size(); k++) {
	//Printf("fill");
        Digit* pDigEvt = &((*digits)[k]);

	int padChX = 0, padChY = 0, module = 0;
    
	o2::hmpid::Digit::pad2Absolute(pDigEvt->getPadID(), &module, &padChX, &padChY);   
       
        if (pDigEvt->mParticlePdg > 10000) {
            hParticlePdgMapAbove[module]->Fill(padChX, padChY);
        } else {
            hParticlePdgMapBelow[module]->Fill(padChX, padChY);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Particle Pdg Map", 800, 600);
    c1->Divide(3, 3); 
    Int_t pos[]= {9,8,6,5,4,2,1};

    for (int i = 0; i < 7; i++) {
        c1->cd(pos[i]);
        hParticlePdgMapBelow[i]->Draw("Colz");
        hParticlePdgMapAbove[i]->Draw("same");
    }
    c1->Show();
    f->Close();
}

void readClusters(char* filename, int nEvent) {
    TFile* f = TFile::Open(filename);
    if (!f) {
        std::cout << "Error: Cannot open file!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)f->Get("o2sim");
    std::vector<Cluster>* clusters = nullptr;
    tree->SetBranchAddress("HMPCluster", &clusters);

    // You would have similar code as readDigits for clusters, but I am leaving it empty for now as you have not provided specifics.
    //...

    f->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <filename> <event_number>" << std::endl;
        return 1;
    }

    char* filename = argv[1];
    int nEvent = std::stoi(argv[2]);

    readDigits(filename, nEvent);
    readClusters(filename, nEvent);

    return 0;
}

