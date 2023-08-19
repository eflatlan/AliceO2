#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)





#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "HMPIDReconstruction/Clusterer.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h> // gRoot

#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fairlogger/Logger.h>
#include <fstream>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"

// C++ header files and libraries
#include <chrono>
#include <ctime>
#include <fstream>
#include <gsl/gsl>
#include <iostream>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <string.h>
#include <ctype.h>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include "randutils.hpp"

#include <queue>


#endif




int range;    // = numTriggers - 1;
int numEvent; // = rng.(0, range)


using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger,
    o2::hmpid::Clusterer;
using std::vector, std::cout, std::cin, std::endl;
using std::this_thread::sleep_for;

using o2::InteractionRecord;

std::priority_queue<int, std::vector<int>, std::greater<int>> topSizes, topSizesInd, topSizesDig;


// 
vector<Digit> oneEventDigits; vector<Cluster> oneEventClusters;
//pair<vector<Digit>, > oneEventCluMapsters
vector<vector<Digit>> nestedVectorDig; vector<vector<Cluster>> nestedVectorClu;

std::vector<o2::hmpid::Digit> mDigitsFromFile,
    *mDigitsFromFilePtr = &mDigitsFromFile;
std::vector<o2::hmpid::Trigger> mTriggersFromFile,
    *mTriggersFromFilePtr = &mTriggersFromFile;
// void SaveFolder(string inpFile);
void changeFont();

vector<double> lastTimes, firstTimes, timeOfEvents;
void sortTimes(vector<double> &timeOfEvents, int nEvents); //

double firstTrg, lastTrg;

TH1F *trigSort = new TH1F("Event Time Histogram", "Event Time Histogram", 50,
                          0., 1000000000.);

TH1F *trigSort2 =
    new TH1F("Instantanoeus Event Frequency Histogram",
             "Instantanoeus Event Frequency Histogram", 50, 0., 30000.);

void sortTriggers(vector<Trigger> &sortedTriggers);
// void sortTriggers(vector<Trigger>& sortedTriggers, TGraph& trigTimeSortStd);
double largestDiff = std::numeric_limits<double>::min();
double largestNegDiff = std::numeric_limits<double>::max(); //

vector<string> dig2Clus(const std::string &fileName, vector<Cluster> &clusters,
                        vector<Trigger> &clusterTriggers,
                        vector<Digit> &digits);

bool mReadFile = false;
std::string mSigmaCutPar;
float mSigmaCut[7] = {4, 4, 4, 4, 4, 4, 4};

std::unique_ptr<TFile> mFile; ///< input file containin the tree
std::unique_ptr<TTree> mTree; ///< input tree

int chargeBelow4[7][5];

std::unique_ptr<Clusterer> mRec; // ef: changed to smart-pointer
int trigTimeCount, trigTimeCount2, trigTimeCount3 = 0;

void initFileIn(const std::string &fileName);

void strToFloatsSplit(std::string s, std::string delimiter, float *res,
                      int maxElem = 7);

struct TriggerTimeInf {
  double timeInNs;
  int triggerIndex;

  TriggerTimeInf(double _timeInNs, int _triggerIndex)
      : timeInNs(_timeInNs), triggerIndex(_triggerIndex) {}
};

vector<TriggerTimeInf> triggerInfoVec;

#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

bool padDig[7][160][144] = {{{true}}};
double chargeAvgByEntries[7][160][144] = {
    {{0.}}}; // chargeAvgByEntries, chargeAvgCount
int chargeAvgCount[7][160][144] = {{{0}}};


//void setPadChannel(bool (&padDigOff)[7][160][144], int chamber, int xLow,
//                   int xHigh, int yLow, int yHigh);

void readDigits();
void fillDigMap(vector<Digit> &digits);

void readClusters(int nEvents = 1, bool leadRun = false) {


  if (leadRun) {
    trigSort2->SetBins(50, 0, 2000.);
  }

  changeFont();
  auto folderName = (gSystem->GetWorkingDirectory());

  auto gwd = gSystem->GetWorkingDirectory();
  std::string path = static_cast<string>(gwd);

  auto runNumber = (gwd.substr(gwd.length() - 9, gwd.length()));

  std::array<std::unique_ptr<TGraph>, 7> trigGraph;

  std::array<std::unique_ptr<TH1F>, 7> digCharge, hMipCharge, digPerEvent,
      digChargeSmall;
  std::array<std::unique_ptr<TH2F>, 7> digMap, digMapAvg, digMapSel, test,
      mapCharge4, digMapEntyAvg, mapEntCount, dig0Charge, oneEventDigMap, oneEventCluMap;

  std::array<std::unique_ptr<TH1F>, 3> triggerTimeFreqHist;
  TH1F trigSortHist;

  std::unique_ptr<TH1F> digPerEventFreq;
  std::unique_ptr<TGraph> trigTime;

  // do not fill high values for selected values
  bool padDigOff[7][160][144] = {{{true}}};
  for (int chamber = 0; chamber < 7; chamber++) {
    for (int x = 0; x < 160; x++) {
      for (int y = 0; y < 144; y++) {
        padDig[chamber][x][y] = true;
      }
    }
  }
  



  vector<Cluster> clusters;
  vector<Trigger> clusterTriggers;
  vector<Digit> digits;
  string fname = "";
  vector<string> fileInfo;
  int numTriggers;

  bool fileFound = false;
  for (const auto &entry : fs::directory_iterator(path)) {
    const auto &pathName = static_cast<string>(entry.path());
    // std::cout << pathName << std::endl;
    if (pathName.size() <= 5)
      continue;

    const auto fileName = pathName.substr(pathName.find_last_of("\\/")+1, pathName.length());
    std::cout << " fileName " << fileName << std::endl;

    if (fileName.size() <= 5)
      continue;


    const auto &fileEnd =
        (fileName.substr(fileName.length() - 5));
    std::cout << " fileEnd " << fileEnd << std::endl;
    if (fileEnd != ".root")
      continue;
    
    
    auto fName = fileName;
    
    for(auto &c : fName)  c = std::tolower(static_cast<unsigned char>(c));

    if(fName.find("hmp") == std::string::npos || fName.find("dig") == std::string::npos)
      continue;

    std::cout << " dig2clus " << pathName << std::endl;


    auto tmp = pathName.substr(0, pathName.find_last_of("\\/"));
    auto folderName = tmp.substr(tmp.length() - 6);
    std::cout << " folderName " << folderName << std::endl;

    std::cout << " folderName.find(\\/) " << folderName.find("\\/") << std::endl;
    std::cout << " std::string::npos " << std::string::npos << std::endl;
    if(folderName.find("\\/") != std::string::npos){
      fname = folderName;
    } else{
      fname = "simulation";
    }

    if (true/*folderName.length() == 6*/) {
      if (true/*std::all_of(folderName.begin(), folderName.end(), ::isdigit)*/) {

        std::cout << " folderName " << folderName << std::endl;
        fileFound = true;
        fileInfo = dig2Clus(pathName, clusters, clusterTriggers, digits);

        // can not sort if only 1 trigger
        if (clusterTriggers.size() > 1) {
          sortTimes(timeOfEvents, nEvents);
        }

        if(fileFound) {
        std::cout << "Found file digit = " << digits.size(); }
        numTriggers = clusterTriggers.size();
        char *fn = strdup(folderName.c_str());
        std::cout << " fname " << fn << std::endl;
        // readDigits();
      } 
    }
    // cout << " Digits file << with {} Clusters, {} Triggers", ,
    // clusters.size(), clusterTriggers.size());
  }
  if (!fileFound) {
    cout << "No fitting file found!";
    return;
    std::exit(0);
  }

  //nestedVectorDig.resize(numTriggers, vector<Digits>);
  //nestedVectorClu.resize(numTriggers, vector<Cluster>);

  //randutils::mt19937_rng rng;
  range = numTriggers - 1;
  numEvent = static_cast<int>(0);//;rng.uniform(0, range);

  float avgDigits = static_cast<float>(1.0f * digits.size() / numTriggers);

  for (int i = 0; i < 7; i++) {


    // ========== MIP Charge Clusters [1D hist 200...2200]=========================
    // canvas = canvas[0
    const char *canStringMip =
        Form("MIP-Charge %i; Charge (ADC channel); Entries/40 ADC", i);
    hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 50, 200., 2200.));
    hMipCharge[i]->SetStats(kTRUE);


    // ========== Digit Charge [1D Histogram 0..2200], Logaritmic x and y-axes=========================
    // canvas[1]
    const char *canStringSize = Form(
        "Logaritmic Digit Charge %i; Charge (ADC channel); Entries/40 ADC", i);
    digCharge[i].reset(new TH1F(canStringSize, canStringSize, 50, 0., 400.));
    digCharge[i]->SetLabelOffset(0.0065, "y");
    digCharge[i]->SetTitleSize(digCharge[i]->GetTitleSize("x") * 1.2, "xy");

    // ========== Digit Charge Small Range [1D Histogram 0..10]=========================
    // canvas = digSmallCanv
    const char *canStringSizes =
        Form("Digit Charge %i;Charge (ADC channel);Entries/40 ADC", i);
    digChargeSmall[i].reset(new TH1F(canStringSizes, canStringSizes, 10, 0., 10.));
    digChargeSmall[i]->SetLabelOffset(0.0065, "y");
    digChargeSmall[i]->SetTitleSize(digChargeSmall[i]->GetTitleSize("x") * 1.2, "xyz");
    digChargeSmall[i]->SetLabelSize(digChargeSmall[i]->GetLabelSize("x") * 1.2, "xyz");


    // ========== Digit MAP of Total Digit Charge [2D Hist 144x160]=========================
    // canvas[2]
    const char *canDigMap = Form("Digit Map %i;x [cm];y [cm]", i);
    digMap[i].reset(new TH2F(canDigMap, canDigMap, 160, 0, 159, 144, 0, 143));
    digMap[i]->SetTitleSize(digMap[i]->GetTitleSize("x") * 1.3, "xyz");
     
    // ========== Digit Map of Charge Below 4 [2D Hist 160x144]=========================
    // canvas digMapLowCan
    const char *mapCharge4Str =
        Form("Chamber %i Digits with Charge < 4 ;x [cm];y [cm]", i);
    mapCharge4[i].reset(
        new TH2F(mapCharge4Str, mapCharge4Str, 160, 0, 159, 144, 0, 143));
    mapCharge4[i]->SetTitleSize(mapCharge4[i]->GetTitleSize("x") * 1.2, "xy");


    // ========== Map of Digits with Charge 0  [2D Hist 160x144] =========================
    // canvas digMap0
    const char *can0Charge = Form("Number of Zero Charge %i; x [cm];y [cm]", i);
    dig0Charge[i].reset(
        new TH2F(can0Charge, can0Charge, 160, 0, 159, 144, 0, 143));
    dig0Charge[i]->SetTitleSize(dig0Charge[i]->GetTitleSize("x") * 1.3, "xy");


    // ========== Map of One Event Digir Hits  [2D Hist 160x144] =========================
    const char *canoneEventDigMap = Form("One Event Digits %i; x [cm];y [cm]", i);
    oneEventDigMap[i].reset(
        new TH2F(canoneEventDigMap, canoneEventDigMap, 160, 0, 159, 144, 0, 143));
    oneEventDigMap[i]->SetTitleSize(oneEventDigMap[i]->GetTitleSize("x") * 1.3, "xy");

    // ========== Map of One Event Cluster Hits  [2D Hist 160x144] =========================
    const char *canOneEventClu = Form("One Event Clusters %i; x [cm];y [cm]", i);
    oneEventCluMap[i].reset(
        new TH2F(canOneEventClu, canOneEventClu, 160, 0, 159, 144, 0, 143));
    oneEventCluMap[i]->SetTitleSize(oneEventCluMap[i]->GetTitleSize("x") * 1.3, "xy");
    
    // ========== Map of Selected Pads by User to Evaluate [2D Hist 160x144]=========================
    // canvas digMapSelCanv
    const char *canDigSel = Form("Selected Pads %i; x [cm];y [cm]", i);
    digMapSel[i].reset(
        new TH2F(canDigSel, canDigSel, 160, 0, 159, 144, 0, 143));
    digMapSel[i]->SetTitleSize(digMapSel[i]->GetTitleSize("x") * 1.3, "xy");


    // ========== Digit MAP of Average Digit Charge Normalized to Number of Events [2D Hist 144x160]=========================
    // canvas canvas[3]
    const char *canDigAvg = Form("Average Charge (Normalized to Total Number "
                                 "of Events) %i;  x [cm];y [cm]",
                                 i);
    digMapAvg[i].reset(
        new TH2F(canDigAvg, canDigAvg, 160, 0, 159, 144, 0, 143));
    digMapAvg[i]->SetTitleSize(digMapAvg[i]->GetTitleSize("x") * 1.3, "xy");


    // ========== Digit Occupancy per Chamber [1D histogram 0.. 0.5%] y [number of Digit-Entries] is logaritmic=========================
    // canvas[4]
    const char *digEvtFreqStr =
        Form("Chamber %i Occupancy ;Occupancy [%];Number of Entries", i);
    digPerEvent[i].reset(new TH1F(digEvtFreqStr, digEvtFreqStr, 500, 0., .5));
    digPerEvent[i]->SetTitleSize(digPerEvent[i]->GetTitleSize("x") * 1.2,
                                 "xyz");
    digPerEvent[i]->SetLabelSize(digPerEvent[i]->GetLabelSize("x") * 1.2, "xy");


    // ================= Average Digit Charge Map by number of Entries per Channel ( [2D histogram 144x160]) ===================
    // canvas digMapEntyAvg
    const char *canEntAvg =
        Form("Average Charge Entry Per Channel %i;  x [cm];y [cm]", i);
    digMapEntyAvg[i].reset(
        new TH2F(canEntAvg, canEntAvg, 160, 0, 159, 144, 0, 143));
    digMapEntyAvg[i]->SetTitleSize(digMapEntyAvg[i]->GetTitleSize("x") * 1.3,
                                   "xyz");


    // ================= Nonzero Digit Entries ( [2D histogram 144x160 mapEntCount]) ===================
    // canvas canvDigitCnt
    const char *canEntCnt =
        Form("Digit-Entires Per Channel %i;  x [cm];y [cm]", i);
    mapEntCount[i].reset(
        new TH2F(canEntCnt, canEntCnt, 160, 0, 159, 144, 0, 143));
    // digMap[i].reset(new TH2F(canDigMap, canDigMap, 160*0.8, 0, 159*0.8,
    // 144*0.8, 0, 160*0.8));
    mapEntCount[i]->SetTitleSize(mapEntCount[i]->GetTitleSize("x") * 1.3,
                                 "xyz");
  }

  for (int chamber = 0; chamber < 7; chamber++) {
    for (int x = 0; x < 160; x++) {
      for (int y = 0; y < 144; y++) {
        if (!padDig[chamber][x][y]) {
          // digMapSel[chamber]->Fill(x, y, 500000.);
          // cout << "False padDigOff " << chamber << " x " << x << " y " << y
          // << endl;
        } else {
          // digMapSel[chamber]->Fill(x, y, 10.);
        }
      }
    }
  }

  for (int chamber = 0; chamber < 7; chamber++) {
    for (int x = 0; x < 160; x++) {
      for (int y = 0; y < 144; y++) {

        // verify that pads are selected off:
        if (digMapSel[chamber]->GetBinContent(x, y) == 50.) {
          // cout << "False padDigOff " << chamber << " x " << x << " y " << y
          // << endl;
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {

    if (i == 0) {
      const char *trigTimeStr = Form("Trigger Time Freq%i", i);
      triggerTimeFreqHist[i].reset(
          new TH1F(trigTimeStr, trigTimeStr, 50, 0, largestDiff));
    } else if (i == 1) {
      const char *trigTimeStr = Form("Instantaneous Event Delta-Time [nS]");
      triggerTimeFreqHist[i].reset(
          new TH1F(trigTimeStr, trigTimeStr, 10, 0, largestDiff));
      triggerTimeFreqHist[i]->SetXTitle("Instantaneous Event Delta-time [nS]");
    } else if (i == 2) {
      const char *trigTimeStr = Form("Instantaneous Event Frequency");
      triggerTimeFreqHist[i].reset(
          new TH1F(trigTimeStr, trigTimeStr, 30000, 0, 30000));
      triggerTimeFreqHist[i]->SetXTitle("Instantaneous Event Frequncy [Hz]");
    }
    triggerTimeFreqHist[i]->SetTitleSize(
        triggerTimeFreqHist[i]->GetTitleSize("x") * 1.2, "xy");
    triggerTimeFreqHist[i]->SetYTitle("Number Of Entries");
  }

  const int nBins = static_cast<int>((lastTrg - firstTrg) * pow(10, -4));
  cout << "Number of Bins " << nBins << endl;
  cout << "1 " << firstTrg << " 2 " << lastTrg << endl;
  trigSort->SetBins(nBins, firstTrg, lastTrg);
  trigSort->SetMaximum(2.);

  const char *trigEvtStr =
      Form("Graph of Instantaneous Event Delta-time; Event Time in LHC nS "
           "(Time Event Occured); Instantaneous Event Delta Time [nS]");
  trigTime.reset(new TGraph);
  trigTime->SetTitle(trigEvtStr);

  std::array<std::unique_ptr<TCanvas>, 8> canvas;
  const char* folder = fname.c_str();
  (canvas[0]).reset(new TCanvas(Form("MIP Cluster-Charge %s", folder),
                                Form("MIP Cluster-Charge %s", folder), 1200,
                                1200));

  (canvas[1]).reset(new TCanvas(Form("Digit Charge Log%s", folder),
                                Form("Digit Charge Log %s", folder), 1200,
                                1200));
  (canvas[1])->SetLogy();
  (canvas[1])->SetLogx();

  (canvas[2]).reset(new TCanvas(Form("Digit-Map %s", folder),
                                Form("Digit-Map %s", folder), 1200, 1200));

  (canvas[3]).reset(new TCanvas(
      Form("Digit-Map Average (Total Events) %s", folder),
      Form("Digit-Map Average (Total Events) %s", folder), 1200, 1200));

  (canvas[4]).reset(new TCanvas(Form("Chamber Occupancy %s", folder),
                                Form("Chamber Occupancy %s", folder), 1200,
                                1200));
  (canvas[4])->SetLogy();

  (canvas[5]).reset(new TCanvas(Form("PlaceHolder %s", folder),
                                Form("PlaceHolder %s", folder), 1200, 1200));
  (canvas[6]).reset(new TCanvas(Form("PlaceHolder2 %s", folder),
                                Form("PlaceHolder2 %s", folder), 1200, 1200));
  (canvas[7]).reset(new TCanvas(Form("PlaceHolder3 %s", folder),
                                Form("PlaceHolder3 %s", folder), 1200, 1200));

  canvas[2]->SetLeftMargin(.1 + canvas[2]->GetLeftMargin());
  canvas[3]->SetLeftMargin(.1 + canvas[3]->GetLeftMargin());
  canvas[4]->SetLeftMargin(.1 + canvas[4]->GetLeftMargin());

  Int_t nTotTriggers = 0;


  std::unique_ptr<Trigger> pTgr;
  std::unique_ptr<Cluster> pClu, pCluEvt, cluster;
  const int digSize = digits.size();
  Printf("digit size = %i", digSize);
  int padChX = 0, padChY = 0, module = 0;

  int trigNum = 0;
  double tprev;
  Trigger trigPrev;

  std::array<float, 7> avgDig;

  const auto rel = 100.0 / (144 * 160);

  for (const auto &trig : mTriggersFromFile) {

    const int numDigPerTrig = trig.getNumberOfObjects();
    const int firstTrig = trig.getFirstEntry();
    const int lastTrig = trig.getLastEntry();

    const auto &orbit = trig.getOrbit();
    const auto &bc = trig.getBc();
    const auto &time = InteractionRecord::bc2ns(bc, orbit);

    const auto &tDif = (trig.getIr()).differenceInBCNS(trigPrev.getIr());

    if (trigNum > 0 && time > pow(10, 6)) {
      const auto &freq = (pow(10, 9)) / tDif;
      triggerTimeFreqHist[2]->Fill(freq);
      triggerTimeFreqHist[1]->Fill(tDif);
      triggerTimeFreqHist[0]->Fill(tDif);
      trigTime->SetPoint(trigNum - 1, static_cast<double>(time), tDif);
      trigSort->Fill(time);
    }

    if (tDif > pow(10, 8)) {
      cout << " Exceeded limit Trigger Number " << trigNum << " Tdiff " << tDif
           << endl;
      cout << "Count " << trigTimeCount3 << endl;
      trigTimeCount3++;
    }

    std::array<int, 7> cntCh = {0, 0, 0, 0, 0, 0, 0};



		oneEventDigits.clear();
    for (int j = firstTrig; j < lastTrig; j++) {
      const auto &dig = digits[j];
      Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);
      cntCh[module]++;
      avgDig[module]++;
			

    oneEventDigits.push_back(dig);
    }
    auto digSize =   lastTrig - firstTrig;

		if(topSizesDig.size() < 10) {
				topSizesDig.push(digSize);
		} else if (digSize > topSizesDig.top()) {
				topSizesDig.pop();
				topSizesDig.push(digSize);
		}


    for (int ch = 0; ch < 7; ch++) {
      digPerEvent[ch]->Fill(rel * cntCh[ch]);
    }

    tprev = time;
    trigPrev = trig;
    trigNum++;
  }

  for (const auto &dig : digits) {
    const auto& digId = dig.getPadID();

    //cout << "digit ID : " << digId << endl;
    Digit::pad2Absolute(digId, &module, &padChX, &padChY);
    const auto &charge = dig.getQ();
    if (padDig[module][padChX][padChY] == true) {
      digCharge[module]->Fill(charge);
      if (charge < 4) {
        chargeBelow4[module][static_cast<int>(charge)] += 1;
        mapCharge4[module]->Fill(padChX, padChY, charge);
      }
      if (charge <= 10) {
        digChargeSmall[module]->Fill(charge);
      }
    }

    digMap[module]->Fill(padChX, padChY, charge);
    digMapAvg[module]->Fill(padChX, padChY, charge / numTriggers);
    if (charge == 0) {
      cout << " charge == 0 " << module << endl;
      dig0Charge[module]->Fill(padChX, padChY, 1.);
    }

    if (charge > 0) {
      chargeAvgByEntries[module][padChX][padChY] += charge;
      chargeAvgCount[module][padChX][padChY]++;

      // digMapEntyAvg[module]->Fill(padChX, padChY, charge);
      mapEntCount[module]->Fill(padChX, padChY, 1.0);
    }

    // digMap[module]->Fill(padChX, padChY, padDigOff[module][padChX][padChY]);
    // padDigits[module][padChX][padChY] += dig.getQ();
  }

  const int clusterSize = clusters.size();
  Printf("clusters size = %i", clusterSize);

  float minCharge = 9999.0f;
  float maxCharge = 0.0f;

  for (auto &clus : clusters) {
    const auto &charge = clus.q();
    const auto &clusSize = clus.size();
    const auto &x = clus.x();
    const auto &y = clus.y();
    const auto &module = clus.ch();
    
    /*
    if (charge == 0) {
      cout << " chargeCluster == 0 | Chamber " << module << endl;
      // digMapSel[module]->Fill(padChX, padChY, 1.);
    } */
     

    if (clusSize >= 3 && clusSize <= 7) {
      hMipCharge[module]->Fill(charge);
    }
  }

  nTotTriggers += clusterTriggers.size();
  const int triggerSize = clusterTriggers.size();
  Printf("trigger size from clusters = %i", triggerSize);

  // auto folderName = fname.c_str();
  std::array<std::unique_ptr<TPaveText>, 7> tpvs;
  std::array<std::unique_ptr<TPaveText>, 7> tpvs2;

  // fileInfo
  const auto f1 = (fileInfo[0]).c_str();
  const auto f2 = (fileInfo[1]).c_str();
  const auto f3 = (fileInfo[2]).c_str();
  const auto f4 = (fileInfo[3]).c_str();

  // const char* runLabel = Form("%i  Duration = %s", fname, f1);
  const char *runLabel = Form("%s", folder);


  // To add textboxes
  vector<const char *> tpvTexts{"MIP Clusters Charge",
                                "Digits-Charge, logx logy",
                                "Digits-Map",
                                "Digits-Map Avg",
                                "Digits Per Event",
                                "Map of Evaluated Areas",
                                "Digits Charge Small Scale"};

  vector<const char *> tpvTexts2{"Event Info",
                                 "Event Information",
                                 "Map of Evaluated Areas",
                                 "Charge Below 4",
                                 "Average Charge Normalized to Channel-Entries",
                                 "Number of Entries in Channel",
                                 "Digits With 0 - charge"};

  int j = 0;
  for (auto &tpv : tpvs2) {
    tpv.reset(new TPaveText(0.05, .05, .9, .9));
    tpv->AddText(Form("%s %s", runLabel, tpvTexts2[j++]));
  }

  for (auto &tpv : tpvs2) {
    tpv->AddText(f1);
    tpv->AddText(f2);
    tpv->AddText(f3);
    tpv->AddText(f4);
  }

  j = 0;
  for (auto &tpv : tpvs) {
    tpv.reset(new TPaveText(0.05, .05, .9, .9));
    tpv->AddText(Form("%s %s", runLabel, tpvTexts[j++]));
  }

  for (auto &tpv : tpvs) {
    tpv->AddText(f1);
    tpv->AddText(f2);
    tpv->AddText(f3);
    tpv->AddText(f4);
  }

  const int tpvSize = static_cast<int>(tpvs.size());
  for (int i = 0; i < tpvSize - 1; i++) {
    canvas[i]->Divide(3, 3);
    canvas[i]->cd(3);
    tpvs[i]->Draw();
  }


  // location of Pads in Canvas to follow HMPID-modules
  static constexpr int posArr[] = {9, 8, 6, 5, 4, 2, 1};



  /// change default font
  changeFont();
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.925);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);


    // ========== Map of One Event Dig-Clu Hits  [2D Hist 160x144] =========================
  const char* oneEvStr = Form("One Event Canvas%s", folder);
  //auto oneEventCanvas = std::make_unique<TCanvas>(oneEvStr, oneEvStr, 1200, 1200);

  std::unique_ptr<TCanvas> oneEventCanvas;// = std::make_unique<TCanvas>
  oneEventCanvas.reset(new TCanvas(oneEvStr, oneEvStr, 1200, 1200));

  oneEventCanvas->Divide(3, 3);
  for(const auto& dig : oneEventDigits){
    const auto& digId = dig.getPadID();
    Digit::pad2Absolute(digId, &module, &padChX, &padChY);
    oneEventDigMap[module]->Fill(padChX, padChY, dig.getQ());
  } 

  Printf("oneEventClusters size = %zu", oneEventClusters.size());
  Printf("oneEventDigits size = %zu", oneEventDigits.size());

  for(auto& clu : oneEventClusters){
    auto chamber = clu.ch();
    auto x = clu.x();
    auto y = clu.y();
    auto charge = clu.q();
    oneEventCluMap[chamber]->Fill(padChX, padChY, charge);
  }

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    TPad *pad = static_cast<TPad *>(oneEventCanvas->cd(pos));
    oneEventDigMap[iCh]->SetStats(kFALSE);
    //oneEventCluMap[iCh]->SetMarkerStyle(3);
    oneEventDigMap[iCh]->Draw("Colz");

    /*oneEventCluMap[iCh]->SetStats(kFALSE);
    oneEventCluMap[iCh]->SetMarkerStyle(3);
    oneEventCluMap[iCh]->DrawCopy("Colz");*/
  }
  oneEventCanvas->SaveAs(Form("One Event Map%s", folder));
  oneEventCanvas->Show();

  // ================= Average Digit Charge Map by number of Entries per Channel (canvDigitAvgEvt [2D histogram 144x160]) ===================
  std::unique_ptr<TCanvas> canvDigitAvgEvt, canvDigitCnt;
  canvDigitAvgEvt.reset(new TCanvas(Form("Event Channel Average %s", folder),
                          Form("Event Channel Average %s", folder), 1200, 2000));
  canvDigitAvgEvt->Divide(3, 3);
  canvDigitAvgEvt->cd(3);
  tpvs2[4]->Draw();

  // ================= Nonzero Digit Entries (canvDigitCnt [2D histogram 144x160 mapEntCount]) ===================
  canvDigitCnt.reset(new TCanvas(Form("Event Channel Count %s", folder),
                          Form("Event Channel Count %s", folder), 1200, 2000));
  canvDigitCnt->Divide(3, 3);
  canvDigitCnt->cd(3);
  tpvs2[5]->Draw();

  int maxCnt[7] = {0};
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    TPad *pad = static_cast<TPad *>(canvDigitAvgEvt->cd(pos));
    for (int x = 0; x < 160; x++) {
      for (int y = 0; y < 144; y++) {
        if (chargeAvgCount[iCh][x][y] > maxCnt[iCh]) {
          maxCnt[iCh] = chargeAvgCount[iCh][x][y];
        }
        auto chAverage =
            chargeAvgByEntries[iCh][x][y] / chargeAvgCount[iCh][x][y];
        if (chargeAvgCount[iCh][x][y] != 0) {
          digMapEntyAvg[iCh]->Fill(x, y, chAverage);
        } else {
          digMapEntyAvg[iCh]->Fill(x, y, 0.);
        }
      }
    }
    cout << "maxCnt " << iCh << maxCnt[iCh] << endl;

    pad->SetLeftMargin(+.035 + pad->GetLeftMargin());
    pad->SetBottomMargin(.085 + pad->GetBottomMargin());
    pad->SetRightMargin(.085 + pad->GetRightMargin());
    pad->SetLogz(1);
    digMapEntyAvg[iCh]->SetStats(kFALSE);
    digMapEntyAvg[iCh]->Draw("Colz");

    TPad *pad2 = static_cast<TPad *>(canvDigitCnt->cd(pos));
    pad2->SetLogz(1);
    pad2->SetLeftMargin(pad->GetLeftMargin());
    pad2->SetBottomMargin(pad->GetBottomMargin());
    pad2->SetRightMargin(pad->GetRightMargin());
    mapEntCount[iCh]->SetStats(kFALSE);
    mapEntCount[iCh]->Draw("Colz");
  }
  canvDigitAvgEvt->SaveAs(Form("AverageByEntries_%s_.png", folder));
  canvDigitCnt->SaveAs(Form("CountEntries_%s_.png", folder));


 
  // ================= Event Time information ===================
  std::unique_ptr<TCanvas> temp1;
  temp1.reset(new TCanvas(Form("Event Information %s", folder),
                          Form("Event Information %s", folder), 1200, 2000));
  temp1->Divide(2, 2);
  temp1->cd(2);
  tpvs2[0]->Draw();

  for (int i = 0; i < 3; i++) {
    TPad *pad2;
    if (i < 1) {
      pad2 = static_cast<TPad *>(temp1->cd(i + 1));
      // largestNegDiff largestDiff
      trigTime->SetMinimum(largestNegDiff);
      trigTime->SetMaximum(largestDiff);
      pad2->SetLeftMargin(.01 + pad2->GetLeftMargin());
      pad2->SetLogy(1);
      // trigTime->SetTitleOffset(trigTime->GetTitleOffset("y")*0.8, "y");
      trigTime->Draw("A*");
    } else {
      pad2 = static_cast<TPad *>(temp1->cd(i + 2));
      pad2->SetLeftMargin(.0275 + pad2->GetLeftMargin());
      pad2->SetBottomMargin(.0375 + pad2->GetBottomMargin());
      pad2->SetRightMargin(.0375 + pad2->GetRightMargin());
      triggerTimeFreqHist[i]->SetTitleOffset(
          triggerTimeFreqHist[i]->GetTitleOffset("y") * 1.2, "xy");

      triggerTimeFreqHist[i]->SetTitleSize(
          triggerTimeFreqHist[i]->GetTitleSize("x") * 0.95, "xy");
      triggerTimeFreqHist[i]->SetLabelSize(
          triggerTimeFreqHist[i]->GetLabelSize("x") * 0.925, "xy");
      triggerTimeFreqHist[i]->Draw();
    }
  }

  gStyle->SetOptStat("eim");
  gStyle->SetStatX(0.95);
  temp1->Show();
  temp1->SaveAs(Form("Event Information_%s_.png", folder));

  /*
   **********************************
   Sorted Triggers
   **********************************
  */

  std::unique_ptr<TCanvas> trigFreqCanv;
  trigFreqCanv.reset(new TCanvas(Form("Trigger Frequency%s", folder),
                          Form("Trigger Frequency%s", folder), 1200, 2000));
  trigFreqCanv->Divide(2, 2);
  // temp2->cd(2);
  // tpvs2[1]->Draw();

  {
    auto pad = static_cast<TPad *>(trigFreqCanv->cd(3));
    pad->SetLeftMargin(.01 + pad->GetLeftMargin());
    // trigTime->SetMinimum(pow(10,12));
    // trigTime->Draw("AC*");
    trigTime->SetMinimum(largestNegDiff);
    trigTime->SetMaximum(largestDiff);
    pad->SetLogy(1);
    trigTime->Draw("A*");
  }

  for (int i = 0; i < 2; i++) {
    TPad *pad = static_cast<TPad *>(trigFreqCanv->cd(1 + i));
    // pad->SetLeftMargin(.0575+pad->GetLeftMargin());
    pad->SetBottomMargin(.0575 + pad->GetBottomMargin());
    // pad->SetRightMargin(.0375+pad->GetRightMargin());
    if (i == 0) {
      trigSort->Draw();
    } else {
      trigSort2->Draw();
    }
  }

  trigFreqCanv->SaveAs(Form("Trigger Frequency Hist and Graph %s.png", folder));


  // ========== Digit MAP of Total Digit Charge [2D Hist 144x160]=========================
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];

    auto pad = static_cast<TPad *>(canvas[2]->cd(pos));

    const auto &pTotalDigs =
        static_cast<float>(100.0f * digMap[iCh]->GetEntries() / digSize);

    pad->SetLogz(1);
    pad->SetLeftMargin(+.035 + pad->GetLeftMargin());
    pad->SetBottomMargin(.085 + pad->GetBottomMargin());
    pad->SetRightMargin(.085 + pad->GetRightMargin());
    digMap[iCh]->SetLabelOffset(digMap[iCh]->GetLabelOffset("y") - 0.0015, "y");
    digMap[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("y") - 0.0035, "xy");
    digMap[iCh]->SetMarkerStyle(3);
    digMap[iCh]->SetStats(kFALSE);

    digMap[iCh]->SetTitle(
        Form("Chamber %i Percentage of total = %02.0f", iCh, pTotalDigs));

    digMap[iCh]->Draw("Colz");
  }

  gStyle->SetStatX(0.85);
  gStyle->SetOptStat("e");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6);

    // ========== Digit MAP of Average Digit Charge Normalized to Number of Events [2D Hist 144x160]=========================
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    auto pad = static_cast<TPad *>(canvas[3]->cd(pos));
    const auto &pTotalDigs =
        static_cast<float>(100.0f * digMapAvg[iCh]->GetEntries() / digSize);
    pad->SetBottomMargin(.085 + pad->GetBottomMargin());
    pad->SetRightMargin(.085 + pad->GetRightMargin());
    pad->SetLeftMargin(.035 + pad->GetLeftMargin());
    digMapAvg[iCh]->SetTitle(
        Form("Chamber Avg %i Percentage of total = %02.0f", iCh, pTotalDigs));
    digMapAvg[iCh]->SetStats(kFALSE);
    pad->SetLogz(1);
    digMapAvg[iCh]->Draw("Colz");
  }

  gStyle->SetOptStat("e");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6);
  gStyle->SetStatX(0.85);


  // ========== MIP Charge Clusters [1D hist 200...2200]=========================
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    auto pad0 = static_cast<TPad *>(canvas[0]->cd(pos));
    pad0->SetLeftMargin(+.0035 + pad0->GetLeftMargin());
    pad0->SetRightMargin(-.005 + pad0->GetRightMargin());
    pad0->SetBottomMargin(.065 + pad0->GetBottomMargin());
    hMipCharge[iCh]->Fit("landau", "I"); // I = fit by integral

    hMipCharge[iCh]->SetTitleOffset(0.8, "y");
    hMipCharge[iCh]->SetTitleSize(hMipCharge[iCh]->GetTitleSize("x") * 1.2,
                                  "xyz");
    hMipCharge[iCh]->SetLabelSize(hMipCharge[iCh]->GetLabelSize("x") * 1.1,
                                  "xyz");
    hMipCharge[iCh]->SetLabelOffset(
        hMipCharge[iCh]->GetLabelOffset("y") * 1.2 * 1.1, "xyz");
    hMipCharge[iCh]->Draw();
  }

  gStyle->SetStatX(0.95);
  // drawMipCharge(hMipCharge)

  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6);
  gStyle->SetOptStat("eimr");
  gStyle->SetLabelOffset(0.00525, "y");


  // ========== Digit Occupancy per Chamber [1D histogram 0.. 0.5%] y [number of Digit-Entries] is logaritmic=========================
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    // ========== Digit MAP =========================
    auto pad5 = static_cast<TPad *>(canvas[4]->cd(pos));
    pad5->SetLeftMargin(+.025 + pad5->GetLeftMargin());
    pad5->SetLogy(1);
    pad5->SetBottomMargin(.1 + pad5->GetBottomMargin());
    pad5->SetRightMargin(-.0025 + pad5->GetRightMargin());
    digPerEvent[iCh]->SetTitleOffset(
    digPerEvent[iCh]->GetTitleOffset("y") + 0.025, "xy");
    digPerEvent[iCh]->Draw();
  }
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6);
  gStyle->SetOptStat("eimr");
  gStyle->SetLabelOffset(0.00525, "y");



    // ========== Map of Digits with Charge 0  [2D Hist 160x144] =========================
  std::unique_ptr<TCanvas> digMap0;
  digMap0.reset(new TCanvas(Form("Digits With 0 - charge%s", folder),
                            Form("Digits With 0 - charge%s", folder), 1200,
                            1200));
  digMap0->Divide(3, 3);

  digMap0->cd(3);
  tpvs2[6]->Draw();
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    auto pad3 = static_cast<TPad *>(digMap0->cd(pos));
    pad3->SetBottomMargin(.0035 + pad3->GetBottomMargin());
    pad3->SetLeftMargin(.085 + pad3->GetLeftMargin());
    pad3->SetRightMargin(.0085 + pad3->GetRightMargin());
    dig0Charge[iCh]->SetTitleOffset(1.3, "y");

    dig0Charge[iCh]->SetLabelOffset(digMap[iCh]->GetLabelOffset("y") - 0.0015,
                                    "y");
    dig0Charge[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("y") - 0.0035,
                                    "y");
    dig0Charge[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("x") - 0.0005,
                                    "x");
    dig0Charge[iCh]->SetMarkerStyle(3);
    dig0Charge[iCh]->SetStats(kFALSE);
    dig0Charge[iCh]->Draw("Colz");
  }

  digMap0->Show();
  digMap0->SaveAs(Form("DigitsWithZeroCharge_%s.png",folder));



  // ========== Map of Selected Pads by User to Evaluate [2D Hist 160x144]=========================
  std::unique_ptr<TCanvas> digMapSelCanv;
  digMapSelCanv.reset(new TCanvas(Form("Pads turned off by user%s", folder),
                                  Form("Pads turned off by user%s", folder),
                                  1200, 1200));
  digMapSelCanv->Divide(3, 3);

  digMapSelCanv->cd(3);
  tpvs2[2]->Draw();
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    auto pad3 = static_cast<TPad *>(digMapSelCanv->cd(pos));
    pad3->SetBottomMargin(.0035 + pad3->GetBottomMargin());
    pad3->SetLeftMargin(.085 + pad3->GetLeftMargin());
    pad3->SetRightMargin(.0085 + pad3->GetRightMargin());
    digMapSel[iCh]->SetLabelOffset(digMapSel[iCh]->GetLabelOffset("y") + 0.0015,
                                   "y");
    digMapSel[iCh]->SetTitleOffset(1.3, "y");
    digMapSel[iCh]->SetMarkerStyle(3);
    digMapSel[iCh]->SetStats(kFALSE);

    digMapSel[iCh]->SetLabelOffset(digMap[iCh]->GetLabelOffset("y") - 0.0015,
                                   "y");
    digMapSel[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("y") - 0.0035,
                                   "y");
    digMapSel[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("x") - 0.0005,
                                   "x");
    digMapSel[iCh]->SetMarkerStyle(3);
    digMapSel[iCh]->SetStats(kFALSE);
    digMapSel[iCh]->Draw("Colz");
  }

  digMapSelCanv->Show();


  // ========== Digit Map of Charge Below 4 [2D Hist 160x144]=========================
  std::unique_ptr<TCanvas> digMapLowCan;
  digMapLowCan.reset(new TCanvas(Form("Charge below 4 %s", folder),
                                 Form("Charge Below 4 %s", folder), 1200, 1200));
  digMapLowCan->Divide(3, 3);
  digMapLowCan->cd(3);
  tpvs2[3]->Draw();
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    auto pad = static_cast<TPad*>(digMapLowCan->cd(pos));

    pad->SetBottomMargin(.085 + pad->GetBottomMargin());
    pad->SetLeftMargin(.025 + pad->GetLeftMargin());
    mapCharge4[iCh]->SetLabelOffset(
    mapCharge4[iCh]->GetLabelOffset("y") + 0.0015, "y");
    mapCharge4[iCh]->SetTitleOffset(.9, "y");
    mapCharge4[iCh]->SetStats(kFALSE);
    mapCharge4[iCh]->SetMarkerStyle(3);
    mapCharge4[iCh]->Draw("Colz");
  }

  digMapLowCan->Show();
  digMapLowCan->SaveAs(Form("Digit Map of Low Charge%s.png", folder));



  // ========== Digit Charge Small Range [1D Histogram 0..10]=========================
  std::unique_ptr<TCanvas> digSmallCanv;
  digSmallCanv.reset(new TCanvas(Form("Digits SmallRange%s", folder),
                      Form("Digits SmallRange %s", folder), 1200, 1200));
  digSmallCanv->Divide(3, 3);
  digSmallCanv->cd(3);
  (tpvs.back())->Draw();

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    auto padDigChargeSmall = static_cast<TPad *>(digSmallCanv->cd(pos));
    padDigChargeSmall->SetBottomMargin(.075 + padDigChargeSmall->GetBottomMargin());
    padDigChargeSmall->SetLeftMargin(.1 + padDigChargeSmall->GetLeftMargin());
    digChargeSmall[iCh]->SetLabelOffset(
    digChargeSmall[iCh]->GetLabelOffset("y") + 0.0015, "y");
    digChargeSmall[iCh]->SetTitleOffset(1.45, "y");
    padDigChargeSmall->SetRightMargin(-.0025 + padDigChargeSmall->GetRightMargin());
    digChargeSmall[iCh]->Draw();
  }


  // ========== Digit Charge [1D Histogram 0..2200], Logaritmic x and y-axes=========================
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    auto padDigCharge = static_cast<TPad *>(canvas[1]->cd(pos));
    (canvas[1])->SetLogy();
    (canvas[1])->SetLogx();
    padDigCharge->SetLogy(1);
    padDigCharge->SetLogx(1);
    padDigCharge->SetBottomMargin(.075 + padDigCharge->GetBottomMargin());
    padDigCharge->SetLeftMargin(.05 + padDigCharge->GetLeftMargin());
    digCharge[iCh]->SetLabelOffset(digCharge[iCh]->GetLabelOffset("y") + 0.0015,
                                   "y");
    digCharge[iCh]->SetTitleOffset(.95, "y");
    digCharge[iCh]->SetTitleOffset(.95, "x");
    digCharge[iCh]->SetLabelSize(digCharge[iCh]->GetLabelSize("x") * 1.05,
                                 "xy");
    padDigCharge->SetRightMargin(-.0025 + padDigCharge->GetRightMargin());

    digCharge[iCh]->Draw();

    for (int charge = 0; charge <= 4; charge++) {
      // cout << " Ch, Charge, num " << iCh << " " << charge << " " <<
      // chargeBelow4[iCh][charge] << endl;
    }
  }

  gStyle->SetStatH(0.2);
  gStyle->SetStatX(0.95);
  gStyle->SetOptStat("eim");
  gStyle->SetLabelOffset(0.008, "y");

  canvas[0]->SaveAs(Form("MIP_Cluster_Charge_%s_.png", folder));
  canvas[1]->SaveAs(Form("Digit_Charge_Log%s_.png", folder));
  canvas[2]->SaveAs(Form("Digit_Map_%s_.png", folder));
  canvas[3]->SaveAs(Form("Digit_Map_Avg_%s_.png", folder));
  canvas[4]->SaveAs(Form("DigitsPerEvent_%s_.png", folder));

  digSmallCanv->SaveAs(Form("Digit_Charge_SmallRange%s_.png", folder));
  digMapLowCan->SaveAs(Form("Digits below 4 %s_.png", folder));
  canvas[0]->Show();
  canvas[1]->Show();
  digSmallCanv->Show();
  canvas[2]->Show();
  canvas[3]->Show();
  canvas[4]->Show();



  sleep_for(500000ms);
}

vector<string> dig2Clus(const std::string &fileName, vector<Cluster> &clusters,
                        vector<Trigger> &clusterTriggers,
                        vector<Digit> &digits) {
  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  uint32_t firstTrigger, lastTrigger = 0;

  mRec.reset(new o2::hmpid::Clusterer()); // ef: changed to smart-pointer

  initFileIn(fileName);

  std::cout << "[HMPID DClusterization - run() ] Enter ...";


  int maxCluSize = 0;
  clusters.clear();
  clusterTriggers.clear();
  cout << "[HMPID DClusterization - run() ] Entries  = " << mTree->GetEntries()
       << endl;
  double durMin, durSec = 0.0;
  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
  } else {
    auto entry = mTree->GetReadEntry() + 1;
    assert(entry < mTree->GetEntries());
    mTree->GetEntry(entry);
    int j = 0;
    int tNum = 0;

    for (const auto &trig : *mTriggersFromFilePtr) {

      LOGP(info, "=========== \n trigger number {} \n", tNum++);

      if (trig.getNumberOfObjects()) {
        gsl::span<const Digit> trigDigits{mDigitsFromFilePtr->data() +
                                              trig.getFirstEntry(),
                                          size_t(trig.getNumberOfObjects())};
        const size_t clStart = clusters.size();
        mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, true); // ef:uncomment

			  /*
				for(int iClu = 0; iClu < clusters.size(); iClu++){
					if(clusters.size() > 0 ) {
						//LOGP(info, "readClusters::dig2Clu : Cluster{} ",iClu);
						auto pClu = &clusters[iClu];
					
						const auto digs = pClu->getDigits();
						
						for(const auto& d : digs) {
							const auto partcilePDG = d->mParticlePdg;
							if(partcilePDG == 211 || partcilePDG == 2212 || partcilePDG == 321) {
								LOGP(info, " particlePDG {} trackID {} eventID {}", 	partcilePDG, d->mTrackID, d->mEventID);
							}
						}
		      }
				}	*/

        clusterTriggers.emplace_back(trig.getIr(), clStart,
                                     clusters.size() - clStart);

				auto cluSize = clusters.size() - clStart;

				if(topSizes.size() < 10) {
					topSizes.push(cluSize);
					topSizesInd.push(tNum);
				} else if (cluSize > topSizes.top()) {
					topSizes.pop();
					topSizes.push(cluSize);

					topSizesInd.pop();
					topSizesInd.push(tNum);
				}



				if(cluSize > maxCluSize) { 
					maxCluSize = cluSize;
					LOGP(info, "new maxCluSize {} -- >  size pads {}",maxCluSize/*, clusters[clStart].getNumDigits()*/, clusters[clStart].size());


				}

				LOGP(info, "Trigger number {}, numClusters {}", tNum, clusters.size() - clStart);			

        if(j == numEvent){
          oneEventClusters = clusters;
          for(const auto& dig : trigDigits)
            oneEventDigits.push_back(dig);
          //oneEventCluMapsters = {trigDigits, clusters};
        }
					
      }
      j++;
    }

    cout << " Received " << mTriggersFromFilePtr->size() << " triggers with "
         << mDigitsFromFilePtr->size()
         << " digits -> clusters = " << clusters.size() << endl;

    digits = *mDigitsFromFilePtr;
    if (digits.size() == 0) {
      digits = mDigitsFromFile;
    }
    mDigitsReceived = mDigitsFromFilePtr->size();
    mClustersReceived = clusters.size();
    mTriggersReceived = mTriggersFromFilePtr->size();
  }

  LOGP(info, "maxCluSize = {}", maxCluSize);
	std::cout << "Top 10 cluster sizes:" << std::endl;
std::vector<int> top10Values;
	while(!topSizes.empty()) {
		top10Values.push_back(topSizes.top());
		topSizes.pop();
	}
	// If you want the values in descending order, you can reverse the vector
	std::reverse(top10Values.begin(), top10Values.end());

	// Print the values
	for(auto val : top10Values) {
		  std::cout << val << std::endl;
	}

	std::vector<int> top10Index;
	while(!topSizesInd.empty()) {
		top10Index.push_back(topSizesInd.top());
		topSizesInd.pop();
	}


	LOGP(info, "with indexes: ");
	for(auto val : top10Index) {
		  std::cout << val << std::endl;
	}


  std::cout << "Top 10 Digit Size" << std::endl;
	std::vector<int> top10ValuesDig;
	while(!topSizesDig.empty()) {
		  top10ValuesDig.push_back(topSizesDig.top());
		  topSizesDig.pop();
	}
	// If you want the values in descending order, you can reverse the vector
	std::reverse(top10ValuesDig.begin(), top10ValuesDig.end());

	// Print the values
	for(auto val : top10ValuesDig) {
		  std::cout << val << std::endl;
	}









  const int numTriggers = static_cast<int>(mTriggersReceived);
  const int numDigits = static_cast<int>(mDigitsReceived);
  const int numClusters = static_cast<int>(mClustersReceived);

  const float digClusRatio = static_cast<float>(1.0f * numDigits / numClusters);
  const float digTrigRatio = static_cast<float>(1.0f * numDigits / numTriggers);
  const auto &ratioInfo =
      Form("Dig/Clus = %.2f Dig/Events= %.0f", digClusRatio, digTrigRatio);

  const auto &digClusInfo =
      Form("Digits %i Clusters %i", numDigits, numClusters);

  if (mDigitsFromFilePtr->size() < 2) {
    const float triggerFrequency =
        static_cast<float>(1.0f * numTriggers / durSec);
    const auto &trigInfo = Form("Events %i, Average Frequency [Hz] = %.2f ",
                                numTriggers, triggerFrequency);
    return {trigInfo, digClusInfo, ratioInfo,
            Form("Not enough Events (%i) for Frequency", numTriggers)};
  }

  // sort triggers by time
  sortTriggers(*mTriggersFromFilePtr);

  Trigger trigFirst = mTriggersFromFilePtr->at(0);
  Trigger trigLast = mTriggersFromFilePtr->back();

  const auto &tDif = (trigLast.getIr()).differenceInBCNS(trigFirst.getIr());

  durSec = static_cast<double>((tDif) / 1000000000.0);
  durMin = static_cast<double>((durSec) / 60.0);
  const auto durInfo = Form("Duration of Events = %.2f min", durMin);
  const float triggerFrequency =
      static_cast<float>(1.0f * numTriggers / durSec);
  const auto &trigInfo =
      Form("Events %i, Frequency [Hz] = %.2f ", numTriggers, triggerFrequency);

  int trigNum = 0;
  Trigger trigPrev;

  const auto &fTrig = mTriggersFromFile[0];
  const auto &lTrig = mTriggersFromFilePtr->back();

  const auto &irFirst = fTrig.getIr();
  const auto &irLast = lTrig.getIr();

  const auto &orbitFirst = fTrig.getOrbit();
  const auto &bcFirst = fTrig.getBc();

  const auto &orbitLast = lTrig.getOrbit();
  const auto &bcLast = lTrig.getBc();

  const auto &nsFirst = InteractionRecord::bc2ns(bcFirst, orbitFirst);
  const auto &nsLast = InteractionRecord::bc2ns(bcLast, orbitLast);

  const auto prevTimeNsTrig = nsFirst;
  int bcPrev;
  unsigned int orbitPrev;

  cout << "First, Ns " << nsFirst << " BC " << bcFirst << " Orbit "
       << orbitFirst << endl;
  cout << "Last , Ns " << nsLast << " BC " << bcLast << " Orbit " << orbitLast
       << endl;

  for (const auto &trig : *mTriggersFromFilePtr) {
    const int numDigPerTrig = trig.getNumberOfObjects();
    const int firstTrig = trig.getFirstEntry();
    const int lastTrig = trig.getLastEntry();

    auto tDelta = (trig.getIr()).differenceInBCNS(trigPrev.getIr());

    const auto &orbit = trig.getOrbit();
    const auto &bc = trig.getBc();
    const auto timeTrigger = InteractionRecord::bc2ns(bc, orbit);

    triggerInfoVec.emplace_back(timeTrigger, trigNum);
    trigPrev = trig;
    if (trigNum > 0) {
      tDelta = InteractionRecord::bc2ns(bc, orbit) -
               InteractionRecord::bc2ns(bcPrev, orbitPrev);
    }
    trigNum++;
    bcPrev = bc;
    orbitPrev = orbit;
  }

  //cout << " largest difference " << largestDiff << endl;
  //cout << " largest negative   " << largestNegDiff << endl;

  return {trigInfo, digClusInfo, ratioInfo, durInfo};
}

void initFileIn(const std::string &filename) {

  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  // Create the TFIle
  mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());

  mTree.reset((TTree *)mFile->Get("o2sim"));
  if (!mTree) {
    mTree.reset((TTree *)mFile->Get("o2hmp"));
  }

  if (!mTree) {

    LOG(warn)
        << "HMPID DigitToClusterSpec::init() : Did not find o2sim tree in "
        << filename.c_str() << endl;
    return;
    std::exit(0);
  }

  if ((mTree->GetBranchStatus("HMPDigit")) == 1) {
    mTree->SetBranchAddress("HMPDigit", &mDigitsFromFilePtr);
  } else if ((mTree->GetBranchStatus("HMPIDDigits")) == 1) {
    mTree->SetBranchAddress("HMPIDDigits", &mDigitsFromFilePtr);
  } else {
    LOG(warn) << "HMPID DigitToClusterSpec::init() : Error in branches!"
              << endl;
    return;
    std::exit(0);
  }

  mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);
  mTree->Print("toponly");
}

void setPadChannel(bool (&padDigOff)[7][160][144], int chamber, int xLow,
                   int xHigh, int yLow, int yHigh) {

  if (xLow > xHigh || yLow > yHigh) {
    cout << "Wrong Comparison" << endl;
    return;
  }

  if (xLow > 159 || xLow < 0 || xHigh > 159 || xHigh < 0) {
    cout << "Wrong Comparison" << endl;
    return;
  }

  if (yLow > 144 || yHigh < 0 || yLow > 144 || yHigh < 0) {
    cout << "Wrong Comparison" << endl;
    return;
  }

  for (int x = xLow; x <= xHigh; x++) {
    for (int y = yLow; y <= yHigh; y++) {
      padDigOff[chamber][x][y] = false;
      padDig[chamber][x][y] = false;
    }
  }
}

void sortTimes(vector<double> &timeOfEvents,
               int nEvents) // lastTimes firstTimes
{

  std::sort(timeOfEvents.begin(), timeOfEvents.end(),
            [](const auto &a, const auto &b) { return (a < b); });

  const int timeLimit = timeOfEvents.size() - nEvents;
  int i = 0;
  for (const auto &time : timeOfEvents) {
    if (i < timeLimit) {
      firstTimes.emplace_back(time);
    } else {
      lastTimes.emplace_back(time);
    }

    i++;
  }

  double avgTimeF, avgTimeL;
  for (const auto &f : firstTimes) {
    // cout << "time : " << f << endl;
    avgTimeF += f;
  }
  avgTimeF = avgTimeF / timeLimit;

  for (const auto &f : lastTimes) {
    cout << "time : " << f << endl;
    avgTimeF += f;
  }
  avgTimeL = avgTimeL / nEvents;
  if (firstTimes.size() > 2) {
    cout << "firstTimes[timeLimit-2]" << firstTimes[timeLimit - 2] << endl;
    cout << " min min2 nax times First : " << firstTimes[0] << " "
         << firstTimes[1] << " " << firstTimes.back() << endl;
  }

  if (lastTimes.size() > 2) {
    cout << " min min2 nax times Last : " << lastTimes[0] << " " << lastTimes[1]
         << " " << lastTimes.back() << endl;
  }
  cout << " avg times First : " << avgTimeF << " Last : " << avgTimeL << endl;
}

void sortTriggers(vector<Trigger> &sortedTriggers) {
  cout << "Sorting triggers length =  " << sortedTriggers.size() << endl;

  std::sort(sortedTriggers.begin(), sortedTriggers.end(),
            [](const auto &a, const auto &b) {
              return (a.getIr()).bc2ns() < (b.getIr()).bc2ns();
            });

  cout << "Iterating through Sorted triggers " << endl;

  int trigNum = 0;
  Trigger trigPrev;

  const auto firstTrig = sortedTriggers[0];
  const auto lastTrig = sortedTriggers.back();
  const auto &tDifTotal =
      (lastTrig.getIr()).differenceInBCNS(firstTrig.getIr());

  cout << "Diff Between first and last Event " << tDifTotal << endl;
  cout << " Avg Frequency = "
       << sortedTriggers.size() / (tDifTotal * pow(10, 9)) << endl;

  firstTrg = (firstTrig.getIr()).bc2ns();
  lastTrg = (lastTrig.getIr()).bc2ns();
  for (const auto &trig : sortedTriggers) {

    const auto &tS = (trig.getIr()).bc2ns();
    const auto &tE = (trigPrev.getIr()).bc2ns();

    const auto &tDif2 = tE - tS;
    const auto &tDif = (trig.getIr()).differenceInBCNS(trigPrev.getIr());

    // cout << "  tDif " << tDif << endl;
    // cout << "  tDif2 " << tDif2 << endl;

    if (trigNum > 0) {
      timeOfEvents.emplace_back(tDif);
      const auto &freq = (pow(10, 9)) / tDif;
      trigSort2->Fill(freq);

      if (tDif > largestDiff) {
        largestDiff = tDif;
      }
      if (tDif < largestNegDiff) {
        largestNegDiff = tDif;
      }
    }
    trigPrev = trig;
    trigNum++;
  }

  cout << "Triggers Sorted" << endl;

  trigSort->SetTitleSize(trigSort->GetTitleSize("x") * 1.3, "xyz");
  trigSort->SetLabelSize(trigSort->GetLabelSize("x") * 1.2, "xyz");
  trigSort->SetXTitle("Time of Event in LHC nS");
  trigSort->SetYTitle("Number of Entries");
  trigSort->SetTitleOffset(trigSort->GetTitleOffset("x") * 1.2, "x");

  trigSort2->SetTitleSize(trigSort2->GetTitleSize("x") * 1.3, "xyz");
  trigSort2->SetLabelSize(trigSort2->GetLabelSize("x") * 1.3, "xyz");
  trigSort2->SetXTitle("Instantaneous Event Frequency [Hz]");
  trigSort2->SetYTitle("Number of Entries");
  trigSort2->SetTitleOffset(trigSort2->GetTitleOffset("x") * 1.2, "x");
}

void strToFloatsSplit(std::string s, std::string delimiter, float *res,
                      int maxElem) {
  int index = 0;
  size_t pos_start = 0;
  size_t pos_end;
  size_t delim_len = delimiter.length();
  std::string token;
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res[index++] = std::stof(token);
    if (index == maxElem) {
      return;
    }
  }
  res[index++] = (std::stof(s.substr(pos_start)));
  return;
}

void changeFont() {
  /*
  std::unique_ptr<TStyle> mStyle;
  mStyle.reset(new TStyle("canvasStyle", "Canvas Root Styles"));
  */
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.925);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.25);
  gStyle->SetStatFontSize(0.065);
  gStyle->SetLegendTextSize(0.08); //
  gStyle->SetTitleSize(.055, "xzy");
  gStyle->SetTitleOffset(.925, "xz"); //.95
  gStyle->SetTitleOffset(1, "y");     // 1.1
  gStyle->SetTitleFontSize(.05);
  // gStyle->SetTitleFont(16, "xz");
  gStyle->SetLabelOffset(0.0065, "y");
  gStyle->SetLabelFont(22, "xyz");
  gStyle->SetLabelSize(.055, "xyz"); //.0525 // verdi av akser
  // mStyle->SetStyle("canvasStyle");
}
