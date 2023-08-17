#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include "populate.cpp"

#include "populate2.cpp" // TODO: change name of class and file here 

#include <math.h>
//#include "ReconE.cpp"
#include "ReconG.cpp"

//namespace ParticleUtils
using namespace o2;

struct Bin {
    float x;
    float y;
};


class ArrAndMap
{


  
  public :
  

  	
  	
  using vecArray2 = std::vector<std::array<double,2>>;
  
  
  vecArray2 arrMaxPionPos, arrMaxKaonPos, arrMaxProtonPos;		
  vecArray2 arrMinPionPos, arrMinKaonPos, arrMinProtonPos;		
  	
  	
  void setMaxArrays(const vecArray2& _arrMaxPionPos, const vecArray2& _arrMaxKaonPos, const vecArray2& _arrMaxProtonPos) {
  	arrMaxPionPos = _arrMaxPionPos;
  	arrMaxKaonPos = _arrMaxKaonPos;
  	arrMaxProtonPos = _arrMaxProtonPos;
  }	
		
  	
  void setMinArrays(const vecArray2& _arrMinPionPos, const vecArray2& _arrMinKaonPos, const vecArray2& _arrMinProtonPos) {
    arrMinPionPos = _arrMinPionPos;
  	arrMinKaonPos = _arrMinKaonPos;
  	arrMinProtonPos = _arrMinProtonPos;
  }	
		
  
	void fillMapFromVec(TH2F* map, const vecArray2& vecArr)
	{
		for(const auto vE : vecArr) {

			//Printf("fillMapFromVec El %.2f %.2f", vE[0], vE[1]);
			map->Fill(vE[0], vE[1]);
		}
	}

  int scale = 10;
  /*
	TH2F* ckovCandMapRange = nullptr;

	TH2F* bgCandMapRange  = nullptr;

	TH2F* ckovCandMapOutRange = nullptr;

	TH2F* bgCandMapOutRange = nullptr;


 TH2F* hSignalAndNoiseMap =  nullptr;

	TH2F*	hSignalMIP =  nullptr;
	TH2F*	hSignalMIPpc =  nullptr;

		TH2F* arrayMaxMin[8];
	TH2F*	hMaxProton =  nullptr;
	TH2F*	hMaxPion =  nullptr;
	TH2F*	hMaxPionMinL =  nullptr;
	TH2F*	hMinPionMaxL =  nullptr;
		TH2F* hMaxKaon =  nullptr;


	TH2F*	hMinProton =  nullptr;
	TH2F*	hMinPion =  nullptr;
	TH2F*	hMinKaon = nullptr; 

  Populate2* populatePtr = nullptr; */ 
		

  void fillCkovCandMapRange(double x, double y) {
   	ckovCandMapRange.emplace_back(std::array<double,2>{x,y});
  }
  
  void fillbgCandMapRange(double x, double y) {
   	bgCandMapRange.emplace_back(std::array<double,2>{x,y});
  }

  void fillckovCandMapOutRange(double x, double y) {
   	ckovCandMapOutRange.emplace_back(std::array<double,2>{x,y});
  }
  
  void fillbgCandMapOutRange(double x, double y) {
   	bgCandMapOutRange.emplace_back(std::array<double,2>{x,y});
  }
  

  
   


	vecArray2 bgCandMapOutRange, ckovCandMapRange, bgCandMapRange, ckovCandMapOutRange;
		
		/*
		~ArrAndMap() {
	hSignalMIP =  nullptr;
		hSignalMIPpc =  nullptr;
		hMaxProton =  nullptr;
		hMaxPion =  nullptr;
		hMaxPionMinL =  nullptr;
		hMinPionMaxL =  nullptr;
		hMaxKaon =  nullptr;


		hMinProton =  nullptr;
		hMinPion =  nullptr;
		hMinKaon = nullptr;

			delete ckovCandMapRange;

			delete bgCandMapRange;

		delete 	ckovCandMapOutRange;

		delete	bgCandMapOutRange;


delete 			hSignalAndNoiseMap;

		delete 	hSignalMIP;
		delete 	hSignalMIPpc;


			delete hMaxProton;
			delete hMaxPion;
			delete hMaxPionMinL;
		delete 	hMinPionMaxL;
			delete hMaxKaon;


			delete hMinProton;
			delete hMinPion;
			delete hMinKaon;

			delete populatePtr;
		
			delete ckovCandMapRange;

			delete bgCandMapRange;

			delete ckovCandMapOutRange;

			delete bgCandMapOutRange;


			delete hSignalAndNoiseMap;

			delete hSignalMIP;
			delete hSignalMIPpc;


			delete hMaxProton;
			delete hMaxPion;
			delete hMaxPionMinL;
			delete hMinPionMaxL;
			delete hMaxKaon;


			delete hMinProton;
			delete hMinPion;
			delete hMinKaon;

			delete populatePtr;
			
		}*/

		int eventCnt;
		
	Populate2* populatePtr = nullptr;
	
	void setPopulatePtr(Populate2* _populatePtr) { 
		populatePtr = _populatePtr; 
	}	
	
	int eventCount = 0;
	void setEventCount(int _eventCount) { 
		eventCount = _eventCount; 
	}	


	void drawTotalMap()
	{

		TH2F* hCkovCandMapRange = nullptr;
		TH2F* hBgCandMapRange  = nullptr;
		TH2F* hCkovCandMapOutRange = nullptr;
		TH2F* hBgCandMapOutRange = nullptr;


		TH2F* hSignalAndNoiseMap =  nullptr;
		TH2F*	hSignalMIP =  nullptr;
		TH2F*	hSignalMIPpc =  nullptr;

		TH2F* arrayMaxMin[8];
		TH2F*	hMaxProton =  nullptr;
		TH2F*	hMaxPion =  nullptr;
		TH2F*	hMaxPionMinL =  nullptr;
		TH2F*	hMinPionMaxL =  nullptr;
		TH2F* hMaxKaon =  nullptr;


		TH2F*	hMinProton =  nullptr;
		TH2F*	hMinPion =  nullptr;
		TH2F*	hMinKaon = nullptr; 

		//Populate2* populatePtr = nullptr;

		hCkovCandMapRange = new TH2F("ckovCandMapRange", "ckovCandMapRange", 160*scale, 0, 159, 144*scale, 0, 143);

		hBgCandMapRange  = new TH2F("bgCandMapRange", "bgCandMapRange", 160*scale, 0, 159, 144*scale, 0, 143);

		hCkovCandMapOutRange = new TH2F("ckovCandMapOutRange", "ckovCandMapOutRange", 160*scale, 0, 159, 144*scale, 0, 143);

		hBgCandMapOutRange = new TH2F("bgCandMapOutRange", "bgCandMapOutRange", 160*scale, 0, 159, 144*scale, 0, 143);

		hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);

		hSignalMIP = new TH2F("hmip ", "hmip ; x [cm]; y [cm]",160*scale*5,0.,159.,144*scale*5,0,143);
		hSignalMIPpc = new TH2F("hmip pc", "hmip pc; x [cm]; y [cm]",160*scale*5,0.,159.,144*scale*5,0,143);


		hMaxProton = new TH2F("maxPoss Ckov Proton", "maxPoss Ckov Proton; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxPion = new TH2F("maxPoss Ckov Pion", "maxPoss Ckov Pion; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxPionMinL = new TH2F("maxPoss Ckov Pion min L", "maxPoss Ckov Pion min L; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMinPionMaxL = new TH2F("minPoss Ckov Pion max L", "minPoss Ckov Pion max L; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxKaon = new TH2F("maxPoss Ckov Kaon", "maxPoss Ckov Kaon; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxProton->SetMarkerColor(kGreen+3);
		hMaxKaon->SetMarkerColor(kRed);

		hMinProton = new TH2F("minPoss Ckov Proton", "minPoss Ckov Proton; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMinPion = new TH2F("minPoss Ckov Pion", "minPoss Ckov Pion; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMinKaon = new TH2F("minPoss Ckov Kaon", "minPoss Ckov Kaon; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);


	  Printf("		fillMapFromVec(mArrAndMap->hMaxPion, arrMaxPion);// map, array");
	  fillMapFromVec(hMaxPion, arrMaxPionPos);// map, array
	  fillMapFromVec(hMaxKaon, arrMaxKaonPos);// map, array
	  fillMapFromVec(hMaxProton, arrMaxProtonPos);// map, array

	  fillMapFromVec(hMinPion, arrMinPionPos);// map, array
	  fillMapFromVec(hMinKaon, arrMinKaonPos);// map, array
	  fillMapFromVec(hMinProton, arrMinProtonPos);
	  
	  

		hMaxProton->SetMarkerColor(kGreen+3);
		hMaxKaon->SetMarkerColor(kRed);
		hCkovCandMapRange->SetMarkerColor(kGreen);
		hCkovCandMapOutRange->SetMarkerColor(kGreen + 2);
		hBgCandMapRange->SetMarkerColor(kRed-1);
		hBgCandMapOutRange->SetMarkerColor(kRed+2);

		hCkovCandMapRange->SetMarkerStyle(3);
		hCkovCandMapOutRange->SetMarkerStyle(2);
		hBgCandMapRange->SetMarkerStyle(3);
		hBgCandMapOutRange->SetMarkerStyle(2);




		
		
	  fillMapFromVec(hCkovCandMapRange, ckovCandMapRange);// map, array
	  fillMapFromVec(hCkovCandMapOutRange, ckovCandMapOutRange);// map, array
	  fillMapFromVec(hBgCandMapRange, bgCandMapRange);// map, array
	  fillMapFromVec(hBgCandMapOutRange, bgCandMapOutRange);

		hMinProton->SetMarkerColor(kGreen+3);
		hMinKaon->SetMarkerColor(kRed);
		hMaxPion->SetMarkerColor(kBlue + 4);    // max ckov max L
		hMinPion->SetMarkerColor(kBlue); 				// min ckov min L
		hMinPionMaxL->SetMarkerColor(kBlue); 		// min ckov max L
		hMaxPionMinL->SetMarkerColor(kBlue + 4);// max ckov min L
		/*
		hMaxPion->SetMarkerStyle(3);
		hMinPionMaxL->SetMarkerStyle(3);*/
		hSignalMIPpc->SetMarkerStyle(3);
		hSignalMIP->SetMarkerStyle(3);
		hSignalMIPpc->SetMarkerStyle(3);
		hSignalMIPpc->SetMarkerColor(kRed);
		hSignalAndNoiseMap->SetMarkerStyle(2);


		
		const auto trkRad = populatePtr->getTrackPos();
		const auto trkPC = populatePtr->getPcImp();
		TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
		TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
		trkRadMap->Fill(trkRad.X(), trkRad.Y());
		trkPCMap->Fill(trkPC.X(), trkPC.Y());


		TCanvas* tcnvRane = new TCanvas(Form("tcnvRane%d", eventCnt), Form("tcnvRane%d", eventCnt), 1600, 800);
		tcnvRane->cd();
		hCkovCandMapRange->Draw();
		hCkovCandMapOutRange->Draw("same");
		hBgCandMapRange->Draw("same");
		hBgCandMapOutRange->Draw("same");


		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");


		trkRadMap->Draw("same");
		trkPCMap->Draw("same");
	}
  /*
	void drawMaxRegions()
	{


		const auto trkPC = populatePtr->getPcImp();
		const auto trkRad = populatePtr->getTrackPos();
		TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
		TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
	  trkRadMap->Fill(trkRad.X(), trkRad.Y());
	  trkPCMap->Fill(trkPC.X(), trkPC.Y());


		trkRadMap->SetMarkerStyle(3);  trkRadMap->SetMarkerColor(kGreen+4);
		trkPCMap->SetMarkerStyle(3);  trkPCMap->SetMarkerColor(kGreen+2);

		TCanvas *thSignalAndNoiseMap = new TCanvas(Form("hSignalAndNoiseMap%d", eventCnt),Form("hSignalAndNoiseMap%d", eventCnt),800,800);  
		thSignalAndNoiseMap->cd();
		hSignalAndNoiseMap->Draw();


		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
		
		hSignalMIP->Draw("same");
		hSignalMIPpc->Draw("same");
		
		trkRadMap->Draw("same");
		trkPCMap->Draw("same");
	}



	void drawTotalMapAndMaxRegions()
	{
		TCanvas* tcnvRane = new TCanvas(Form("TotalMap%d", eventCnt), Form("TotalMap%d", eventCnt), 1600, 800);

				const auto trkPC = populatePtr->getPcImp();
		const auto trkRad = populatePtr->getTrackPos();
		TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
		TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
	  trkRadMap->Fill(trkRad.X(), trkRad.Y());
	  trkPCMap->Fill(trkPC.X(), trkPC.Y());


		trkRadMap->SetMarkerStyle(3);  trkRadMap->SetMarkerColor(kGreen+4);
		trkPCMap->SetMarkerStyle(3);  trkPCMap->SetMarkerColor(kGreen+2);

		tcnvRane->cd();
		ckovCandMapRange->Draw();
		ckovCandMapOutRange->Draw("same");
		bgCandMapRange->Draw("same");
		bgCandMapOutRange->Draw("same");

		hSignalAndNoiseMap->Draw();


		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
		
		hSignalMIP->Draw("same");
		hSignalMIPpc->Draw("same");
				trkRadMap->Draw("same");
		trkPCMap->Draw("same");
	} */ 

};

class CkovTools {


private: 

  
	// for storing candidates...
	struct Candidate {
		double x, y = 0.;
		double R = 0.;
		double phiL = 0., phi = 0.;
		bool isCandidate = false;
	};


  bool print = false; // ef: TODO : later pass this in ctor 

  Populate2* populate2Ptr = nullptr;// =  
	Populate2* populatePtr = nullptr;// 

// using array = std::array;
using vecArray4 = std::vector<std::array<double,4>>;
using vecArray3 = std::vector<std::array<double,3>>;
using vecArray2 = std::vector<std::array<double,2>>;

using vecArrayPair = std::vector<std::pair<double,double>>;


// using arrArray3 = std::vector<std::array<double,3>>;


using segType = std::vector<std::pair<double, double>>;
segType segPionLocal = {{}};

bool kaonStatus = true, pionStatus = true, protonStatus = true;
//TLine* tlinePion;

static constexpr double PI = M_PI;
static constexpr double halfPI = M_PI/2;
static constexpr double twoPI = M_PI*2;


static constexpr double stdDevPion = 0.008; 
static constexpr double stdDevKaon = 0.008; 
static constexpr double stdDevProton = 0.008;
static constexpr float tGap = 8;
static constexpr float  rW = 1.5; // was 1?
static constexpr float  qW = 0.5;
static constexpr float lMax = 1;


static constexpr float CH4GapWidth = 8;
static constexpr float  RadiatorWidth = 1.5; // was 1?
static constexpr float  QuartzWindowWidth = 0.5;


const float L = rW/2;


  // L value for reconstruction
  static constexpr float  EmissionLenght = RadiatorWidth/2;

  float thetaP, phiP, xPC, yPC, xRad, yRad; 
 float nF, nQ, nG;  
 std::array<float, 3> ckovHyps;
 std::vector<std::pair<double, double>> photons;


 // instantiate ckovhyp boundaries to "invalid values" 
 // --> later set them if they exceed momentum threshold
 double ckovPionMin = 999, ckovPionMax = 0, ckovKaonMin = 999, ckovKaonMax = 0, ckovProtonMin = 999,ckovProtonMax = 0;

 // double mRMax, mL1Max, mL2Max, mRMin, mL1Min, mL2Min;

 double momentum, mass;

 float cosThetaP, sinThetaP, tanThetaP;
 float cosPhiP, sinPhiP, tanPhiP;
 float trackCkov;
 double xMipLocal, yMipLocal;
 float phiLCurrent, etaCCurrent;


TVector2 trkPos;
TVector3 trkDir; 

 //std::unique_ptr<ArrAndMap> mArrAndMap;// = nullptr;// new ArrAndMap(eventCnt);

 int eventCnt;
 TVector2 trkPC;

public:
typedef std::vector<std::pair<double, double>> MapType;

void addPhoton(double phiL, double etaC)
{
  photons.push_back(std::make_pair(phiL, etaC));
}


// set new values of phiL and etaC for new photon
void setPhoton(double phiL, double etaC)
{
  phiLCurrent = phiL;
  etaCCurrent = etaC;
}


//   double pc[5] = {xRad,yRad,L,thetaP, phiP};
double refIndexes[3] = {nF, nQ, nG};


~CkovTools() {
}


// ef: let trackCkov be empty aon; TODO: add this based on pdg!


CkovTools (double radParams[7], double refIndexes[3], double MIP[3],
           std::array<float, 3> ckovHyps, float trackCkov, int eventCnt)
  : 
    ckovHyps(ckovHyps), , trackCkov(trackCkov), eventCnt(eventCnt) {  
    
    
  // double radParams[6] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum};
  
  xMip = MIP[0], yMip = MIP[1], qMip = MIP[2]; 

  xRad= radParams[0];
  yRad= radParams[1];
  L = radParams[2]; 
  thetaP = radParams[3];
  phiP = radParams[4];
  momentum = radParams[5];	 
  mass = radParams[6]; // ef; let this be empty aon! TODO: get this from pdg

  nF = refIndexes[0];

	// set tehse to be constant?
  nQ = 1.5787; // ??? 1.5787;// TODO: check this !
  nG = 1.005;
 
  Printf(" CkovTools momentum = %.2f, refFreon = %.2f; ckovHyps : %.2f %.2f %.2f", momentum, nF, ckovHyps[0], ckovHyps[1], ckovHyps[2]);

	if(TMath::IsNaN(ckovHyps[0])){
 	  Printf("Pion CkovHyps is Nan!");
	  setPionStatus(false);
	  // setIsNan()?
	} else {
		setPionStatus(true);
 	  Printf("Pion CkovHyps %.2f", ckovHyps[0]);
  }

  if(TMath::IsNaN(ckovHyps[1])){
	  Printf("Kaon CkovHyps is Nan!");
   	setKaonStatus(false);
	  // setIsNan()?
	} 
	else { 
	  	//setMaxRadius();
  }

  if(TMath::IsNaN(ckovHyps[2])){
  	  Printf("Proton CkovHyps is Nan!");
	  	setProtonStatus(false);
	} else {
		  //setMaxRadius();
	}

	if(getPionStatus()){
	  ckovPionMin = ckovHyps[0] - 4 * stdDevPion;
	  ckovPionMax = ckovHyps[0] + 4 * stdDevPion;
 	  Printf("init CkovTools constructor : getPionStatus() true ! minPion %.2f, maxPion %.2f ", ckovPionMin, ckovPionMax);
  }	else {
 	  Printf("init CkovTools constructor : getPionStatus() was false !");
  }
  
  if(getKaonStatus()){
	  ckovKaonMin = ckovHyps[1] - 4 * stdDevKaon;
  	ckovKaonMax = ckovHyps[1] + 4 * stdDevKaon;
  } if(getProtonStatus()){
		ckovProtonMin = ckovHyps[2] - 4 * stdDevProton;
		ckovProtonMax = ckovHyps[2] + 4 * stdDevProton;
	}

cosThetaP = TMath::Cos(thetaP);
sinThetaP = TMath::Sin(thetaP);
tanThetaP = TMath::Tan(thetaP);

cosPhiP = TMath::Cos(phiP);
sinPhiP = TMath::Sin(phiP);


TRotation rotZ; rotZ.RotateZ(phiP);
TRotation rotY; rotY.RotateY(thetaP);

TVector3 ip(0,0,rW-L+tGap+qW);
TVector3 op; op = rotZ*rotY*ip;

xMipLocal = // tanThetaP*cosPhiP*(rW-L + tGap + qW);
yMipLocal = // tanThetaP*sinPhiP*(rW-L + tGap + qW);

auto dX = tanThetaP*cosPhiP*(rW-L + tGap + qW);
auto dY = tanThetaP*sinPhiP*(rW-L + tGap + qW);


xPC = dX + xRad; // bruke PC eller MIP?
yPC = dY + yRad;


// fyll : 


// trkPC.Set(xPC, yPC); // TODO :check it equals populate.getPcImp();
trkPC.Set(xPC, yPC); // TODO :check it equals populate.getPcImp();

 
trkPos.Set(xRad, yRad);
trkDir; 
trkDir.SetMagThetaPhi(1, thetaP, phiP);

populatePtr = new Populate2(trkPos, trkDir, nF, L);
populatePtr->setPcImp(trkPC);

//mArrAndMap = new ArrAndMap(eventCnt, populatePtr);
//mArrAndMap = std::make_unique<ArrAndMap>(eventCnt, populatePtr);



/*
(mArrAndMap->hSignalMIP)->Fill(xRad, yRad);
(mArrAndMap->hSignalMIPpc)->Fill(xPC, yPC);
*/ 

populate2Ptr = new Populate2(trkPos, trkDir, nF, L);

populate2Ptr->setPcImp(trkPC);

Populate populate(trkPos, trkDir, nF, L);

populate.setPcImp(trkPC);


Printf("init Ckovtools \n MIP Root : %f %f %f \n MIP local %f %f",op.Px(),op.Py(),op.Pz(),xMipLocal,yMipLocal);
      // constructor body goes here, if needed

Printf("Ckovtools :: trkPC %.2f %.2f",trkPC.X(), trkPC.Y());
      // constructor body goes here, if needed

/*
mRMax = getR_Lmax(ckovPionMax, halfPI);
mL2Max = getR_Lmax(ckovPionMax, 0);
mL1Max = getR_Lmax(ckovPionMax, PI);


// needs to be changed if to be used! now uses L = lMax
mRMin = getR_Lmax(ckovProtonMin, halfPI);
      mL1Min = getR_Lmax(ckovProtonMin, PI);
      mL2Min = getR_Lmax(ckovProtonMin, 0);*/ 
		
			/*for(const auto& c : ckovHyps) {
				auto R = getR_Lmax(c, halfPI);
				auto l2 = getR_Lmax(c, 0);
				auto l1 = getR_Lmax(c, PI);
				Printf(" CkovTools : CkovHyp %f, R %f, l1 %f l2 %f", c, R, l2, l1);
      } */  
}





void setPionStatus(bool status)
{ 
  pionStatus = status;
}

bool getPionStatus() const
{
  return pionStatus;
}

void setKaonStatus(bool status)
{ 
  kaonStatus = status;
}

bool getKaonStatus() const
{
  return kaonStatus;
}

void setProtonStatus(bool status)
{ 
  protonStatus = status;
}

bool getProtonStatus() const
{
  return protonStatus;
}


double getMinCkovKaon()
{
  return ckovKaonMin;
}

double getMaxCkovKaon()
{
  return ckovKaonMax;
}

double getMinCkovPion()
{
  return ckovPionMin;
}

double getMaxCkovPion()
{
  return ckovPionMax;
}

double getMinCkovProton()
{
  return ckovProtonMin;
}

double getMaxCkovProton()
{
  return ckovProtonMax;
}


double getPhiP()
{
  return phiP;
}

double getThetaP()
{
  return thetaP;
}

// TODO: add nQ
// only consider photons in the correct range:


// TODO: add nQ
// only consider photons in the correct range:


// x, y, etaC


// std::array<int, 4> arrayInfo;

	// ;

//std::vector<std::pair<double, double>> segment(std::vector<std::array<double, 3>>& cherenkovPhotons, vecArray2& pionCands, vecArray2& kaonCands, vecArray2& protonCands, std::array<int, 4>& arrayInfo, std::vector<ParticleUtils::Candidate2>& candCombined, MapType& bins) { 
std::vector<std::pair<double, double>> segment(const std::vector<o2::hmpid::Cluster>& clusterTrack, std::array<int, 4>& arrayInfo, int trackIndex, const std::vector<float>& mipCharges, float mipX, float mipY, const int mcmcTrackPdg)
{ 


  vecArray2 pionCandidates, kaonCandidates, protonCandidates;
  
  
  ArrAndMap mArrAndMap;// new ArrAndMap(eventCnt);
  
  Printf("ckovTools enter  ckovTools.segment"); 	 
  //const auto infString = Form("localRef #Theta_{p}  = %.4f #Phi_{p} = %.4f L = %.2f \n #Theta_{C} = %.4f maxCkov = %.4f ; x [cm]; y [cm]", thetaP,phiP, L,trackCkov,ckovPionMax); 

  //TH2F *localRefMIP = new TH2F("localRefMIP ", infString,800,-40.,-40.,900,-40.,40.);

  
  /*
  TVector2 trkPos(xRad, yRad);
  TVector3 trkDir; 
  trkDir.SetMagThetaPhi(1, thetaP, phiP);
  Populate* populatePtr = new Populate(trkPos, trkDir, nF);
  Populate populate(trkPos, trkDir, nF);*/


  // void setArrayMin(Populate* populate, double etaTRS, vecArray3 inPutVector) 


  const size_t kN = 400;

	//array<array<double, 3> ,kN> arrMaxPion; 

  

	// disse kan brukes istedet for maxKaonVec...?
	vecArray4 arrMaxPion, arrMinPion, arrMinProton, arrMaxProton, arrMinKaon, arrMaxKaon;

	vecArray2 arrMaxPionPos, arrMinPionPos, arrMinProtonPos, arrMaxProtonPos, arrMinKaonPos, arrMaxKaonPos;

	// do not reserve size? NO prob just dont resize :) 
	arrMaxPion.reserve(kN);
	arrMinPion.reserve(kN);
	//arrMaxPion.resize(kN);
	//arrMinPion.resize(kN);

	arrMaxPionPos.reserve(kN);
	arrMinPionPos.reserve(kN);
	//arrMaxPionPos.resize(kN);
	//arrMinPionPos.resize(kN);

	// arrMaxPion.resize(kN);
		
	// check if candidate can be proton (i.e., that it exceeds momentum threshold)
	if(getProtonStatus()) {
		arrMaxProton.reserve(kN);
		arrMinProton.reserve(kN);
		//arrMaxProton.resize(kN);
		//arrMinProton.resize(kN);

		arrMinProtonPos.reserve(kN);
		arrMaxProtonPos.reserve(kN);
		//arrMinProtonPos.resize(kN);
		//arrMaxProtonPos.resize(kN);
				
		// also later add max, i forste runde sjekk at ckov i {minProton, maxPion}
		Printf("calling setArrayMin w getMinCkovProton() = %.2f", getMinCkovProton());
  	setArrayMin(getMinCkovProton(), arrMinProton, arrMinProtonPos, kN);

		Printf("calling setArrayMax w getMaxCkovProton() = %.2f", getMaxCkovProton());
  	setArrayMax(getMaxCkovProton(), arrMaxProton, arrMaxProtonPos, kN);
	}
	

	// check if candidate can be Kaon (i.e., that it exceeds momentum threshold)
	if(getKaonStatus()) {
		arrMaxKaon.reserve(kN);
		arrMinKaon.reserve(kN);				

		//arrMaxKaon.resize(kN);
		//arrMinKaon.resize(kN);	

		arrMaxKaonPos.reserve(kN);
		arrMinKaonPos.reserve(kN);				

		//arrMaxKaonPos.resize(kN);
		//arrMinKaonPos.resize(kN);

		Printf("calling setArrayMin w getMinCkovKaon() = %.2f", getMinCkovKaon());
  	setArrayMin(getMinCkovKaon(), arrMinKaon, arrMinKaonPos, kN);

		Printf("calling setArrayMax w getMaxCkovKaon() = %.2f", getMaxCkovKaon());
  	setArrayMax(getMaxCkovKaon(), arrMaxKaon, arrMaxKaonPos, kN);
	}

	


	Printf("BF : Length of elem vectors : arrMaxPion %zu", arrMaxPion.size());	
  Printf("calling setArrayMax w getMaxCkovPion() = %.2f", getMaxCkovPion());
  setArrayMax(getMaxCkovPion(), arrMaxPion, arrMaxPionPos, kN);
  
  for(const auto& ip : arrMaxPion) {
		const auto& phiL_ = ip[0];
		const auto& phiR_ = ip[1]; 
		const auto& r_ = ip[2];  
		//Printf("setArrayMax() --> checking arrMaxPionPos | : phiL %.2f, phiR %.2f, r %.2f", phiL_, phiR_, r_);
		if(r_ == 0 ) {throw std::invalid_argument("r====??????;");}
  } 
  
	Printf("AFTER : Length of elem vectors : arrMaxPion %zu", arrMaxPion.size());	



	//Printf(" BF : Length of elem vectors : arrMinPion %zu", arrMinPion.size());
  //Printf("calling setArrayMin w getMinCkovPion() = %.2f", getMinCkovPion());
  setArrayMin(getMinCkovPion(), arrMinPion, arrMinPionPos, kN);
	//Printf("AFTER : Length of elem vectors : arrMinPion %zu", arrMinPion.size());
	// fill all the values in the maps 
	
	//Printf("Length of elem vectors : arrMinPion %zu", arrMinPion.size());
	//Printf("Length of elem vectors : arrMaxPion %zu", arrMaxPion.size());
	
	if(true) {

		/*Printf("		fillMapFromVec(mArrAndMap->hMaxPion, arrMaxPion);// map, array");
		fillMapFromVec(mArrAndMap->hMaxPion, arrMaxPionPos);// map, array
		fillMapFromVec(mArrAndMap->hMaxKaon, arrMaxKaonPos);// map, array
		fillMapFromVec(mArrAndMap->hMaxProton, arrMaxProtonPos);// map, array

		fillMapFromVec(mArrAndMap->hMinPion, arrMinPionPos);// map, array
		fillMapFromVec(mArrAndMap->hMinKaon, arrMinKaonPos);// map, array
		fillMapFromVec(mArrAndMap->hMinProton, arrMinProtonPos);// map, array
		*/ 
	} 
  



  // get track impacrt point at PC


	const int scale = 2;


  Printf("ckovTools segment : exit const auto& p : segPionLocal"); 	 

  // TODO: ckovHyps : get std-dev for Theta_ckov of pion kaon and proton from the values theta_i

 // initialize recon with track input params 
 

    // TODO: change this, xMipLocal just placeholder
    // not sure if xPC simply is obtained like this
    double xPC = xMipLocal, yPC = yMipLocal;
    local2GlobalRef(xPC, yPC);


    Printf("dX %f dY %f ", xMipLocal, yMipLocal);

    int numPhotons= 0;


  const auto area = 144*156; 
 
 
  MapType filledBins;
 
  std::vector<std::pair<double, double>> photonCandidates;


 
    double xML = xMipLocal, yML = yMipLocal;

    //int iPhotCount = 0;





	//const size_t numPhotonsTemp = ckovAndBackground.size();
	//std::vector<ParticleUtils::Candidate2> candidatesCombined;
	//candidatesCombined.reserve(numPhotonsTemp);

	
	//Printf("length ckovPhotons %zu length background %zu length total %zu",cherenkovPhotons.size(), backGroundPhotons.size(), ckovAndBackground.size());
	
	// gjør dette mer elegant snre
	
  // ckovAndBackground


  Printf(" enter const auto& photons : ckovAndBackground");


  
  //for(const auto& photons : ckovAndBackground) {
  for(const auto& photons : clusterTrack) 
  {

    const auto& dist = (photons.x() - mipX)*(photons.x() - mipX) + (photons.y() - mipY)*(photons.y() - mipY);
    

    // match track PDG w MIP cluster : 
    if(photons.x() == mipX && photons.y() == mipY && photons.getQ() ==  mipCharge) { // current cluster is mip
      // check if PDG code matches track's PDG-code
      // : photons has field digits w pdg-code; get this and match w track: 

      std::vector<Topology> topologyVector = photons.getClusterTopology();

      if(topologyVector.size() != 0) {
        const auto photonPDG = topologyVector[0].pdg; // check also other indexes?

        if(photonPDG != mcTrackPdg) {
          Printf("photonPDG != mcTrackPdg"); 
          
        }
        else {
          Printf("photonPDG matched mcTrackPdg!"); 

          // the current photon is the MIP corresponding to the current track
          // good! set the truth or something? 
        }
      }

      if()
     
    }


    if(dist < 2) // add small radius where we dont add; this also to not include the MIP for hte current track;
      continue;  // also to not add candidates from other MIPs


    // this means 
    bool skip = false;
    for(const auto& mipCharge : mipCharges ) { 
       if(photons.getQ() ==  mipCharge) { 
        skip = true; 
       } 
    }

    if(skip) 
      continue; // this photon charge is a MIP--> dont consider as candidate 


    // this means current entry should be a MIP, this will not be a candidate


    // ef : TODO get this value from calibration
    //if(photons.getQ() > mipCut) {continue; } // this means current entry should be a MIP, this will not be a candidate

		iPhotCount++;


  	Printf("%d", iPhotCount);
    //const auto& x = photons[0], y = photons[1];
    const auto& x = photons.x(), y = photons.y();
    


    const auto& etaC = photons[2];


    double xL = x, yL = y;
    double xG = x, yG = y;
    /*localRefUnrot->Fill(x,y);
    // transform to phiRing ref-system
    local2PhiRing(xL, yL, xMipLocal, yMipLocal);
    localRef->Fill(xL, yL);*/ 


		/*
    Printf("ckovtools cherenkov photons x > xMaxPhiPi && x < xMaxPhi0 && y > yMaxPhiPi && y < yMaxPhi0");
    Printf("ckovtools cherenkov photons x  %f > xMaxPhiPi %f && x %f < xMaxPhi0 %f && y %f > yMaxPhiPi  %f && y %f < yMaxPhi0 %f",  x, xMaxPhiPi, x , xMaxPhi0 , y , yMaxPhiPi , y , yMaxPhi0);*/
    
	
    //Printf("\nckovtools cherenkov photons x  %f > -mL1Max %f && x %f < mL2Max %f && y %f > -mRMax  %f && y %f < mRMax %f\n",  x, -mL1Max, x , mL2Max , y , -mRMax , y , mRMax);

    double thetaCer, phiCer;
    //local2GlobalRef(xG, yG);
    // double cluX, double cluY, double& thetaCer, double& phiCer
    

    //Printf("CkovTools segment thetaCer %f phiCer %f", thetaCer, phiCer);


    auto xAbs = TMath::Abs(x);
    auto yAbs = TMath::Abs(y);

    //Printf("\nckovtools cherenkov photons xAbs  %f > mL1Max %f && x %f < mL2Max %f && yAbs %f > mRMax  %f && y %f < mRMax %f\n",  x, mL1Max, xAbs , mL2Max , yAbs , mRMax , y , mRMax);




    // not really helpful? 
    if(true){
        //if(x > -mL1Max && x < mL2Max && y > -mRMax  && y < mRMax){

        //double thetaCer, phiCer;
        //reconG.findPhotCkov(xG, yG, thetaCer, phiCer);	
        //auto ckov = thetaCer;


        bool withinRange = true; 



        const TVector2 posPhoton(x, y);

        // find reference  radius to be compared to segmetation radiuses 
        const auto rPhoton = (posPhoton - trkPC).Mod();


        // skal denne vaere trkRAD??
        const auto phiPhoton = (posPhoton - trkPC).Phi(); // ef : skal denne vaere rad?

        const auto pc = populatePtr->getPcImp();

        
        
        /*
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        */ 
    
        // iteration phiL approach init
        // checking if inside outer pionRadius (Lmin, etaMax)

        
        bool isPhotonProtonCand = false, isPhotonKaonCand = false, isPhotonPionCand = false;

        bool isMaxProtonOk = false, isMinProtonOk = false, isMaxKaonOk = false, isMinKaonOk = false, isMaxPionOk = false, isMinPionOk = false;


        Printf("\n\n ============================================================ \n");		
        if(getProtonStatus() and getKaonStatus() and getPionStatus()) {
            
            Printf("Pion%d can be Pion Kaon and Proton from p-hyp \n", iPhotCount);	
            // verify thatrPhoton > rMinProton(@ phiEstimated = phiPhoton)

            Printf("\n\npopulate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, above, arrMinProton, getMinCkovProton());");
            isMinProtonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinProton, getMinCkovProton(), "Proton");

            // check if rPhoton > rMax(@ phiEstimated = phiPhoton)
            if(isMinProtonOk) {
                Printf("isMinProtonOk = true");
                Printf("\n\npopulate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion());");

                isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");
                
                // this means rProtonMin < rPhoton < rMaxPion
                if(isMaxPionOk) {
                    
                    
                    // this means rProtonMax > rPhoton; then also the rPh < rPionMax and rPh < rKaonMax

                    isMaxProtonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxProton, getMaxCkovProton(), "Proton");
                    Printf("===================================="); 
                    if(isMaxProtonOk) {
                        // we have proton-candiate
                        isPhotonProtonCand = true;

                        protonCandidates.push_back(std::array<double,2>{x,y});
                        Printf("Photon%d is a Proton Candidate", iPhotCount); 
                        // isMaxPionOk = true;
                        isMaxKaonOk = true;
                    }	else {
                        Printf("Photon%d not a Proton Candidate", iPhotCount); 
                    }
                    Printf("====================================\n"); 
                    Printf("===================================="); 
                    // this means rPionMin < rPhoton --> and also rKaonMin < rPhoton



                    isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion");
                    if(isMinPionOk) {
                        // we have pion-candiate
                        isPhotonPionCand = true;
                        Printf("Photon%d is a Pion Candidate", iPhotCount); 
                                            
                        pionCandidates.push_back(std::array<double,2>{x,y});
                        
                        isMinKaonOk = true;
                    }	else {
                        Printf("Photon%d not a Pion Candidate", iPhotCount); 
                    }			
                        Printf("===================================="); 				
                    
                    // if rPhoton < rPhotonMin, check if rPhoton > rKaonMin
                    if(!isMinKaonOk) {
                        // isMinKaonOk = populate2Ptr->checkCond()
                    isMinKaonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinKaon, getMinCkovKaon(), "Kaon");
                    }
                    
                    Printf("\n===================================="); 
                    // only check isMaxKaon if isMinKaonOk
                    if(isMinKaonOk) {
                        if(!isMaxKaonOk) {
                            isMaxKaonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxKaon, getMaxCkovKaon(), "Kaon");
                            // isMaxKaonOk = populate2Ptr->checkCond()
                        } 

                        if(isMaxKaonOk and isMinKaonOk) {
                            Printf("Photon%d is a Kaon Candidate", iPhotCount); 
                            // we have kaon cand
                            isPhotonKaonCand = true;
                            kaonCandidates.push_back(std::array<double,2>{x,y});
                        } else {
                            Printf("Photon%d not a Kaon Candidate", iPhotCount); 
                        }
                    } else {
                        Printf("Photon%d not a Kaon Candidate", iPhotCount); 
                    }
                    Printf("====================================\n"); 
                                
                }  // end if isMaxPionOk == true
                            
                // else isMaxPionOk == false
                else {
                    Printf("===============================================================");	
                    Printf("\n Photon%d not a Hadron Candidate :\n rPhoton %.2f > rPionMax",iPhotCount,  rPhoton);
                    Printf("===============================================================");	
                }
                
            }	// end if isMinProtonOk == true
            
            // else isMinProtonOk == false
            // this means  rPhoton < rProtonMin : 
            else {
                Printf("===============================================================");	
                Printf("\n Photon%d not a Hadron Candidate :\n rPhoton %.2f < rProtonMin",iPhotCount, rPhoton);
                Printf("===============================================================");	
            }

        
        } // end if getProtonStatus() and getKaonStatus and getPionStatus

        // if Proton not is possible because momentum-threshold is not exceeded:
        else if(!getProtonStatus() and getKaonStatus() and getPionStatus()) {
            Printf("\n====================================================="); 
            Printf("Photon%d can be Pion and Kaon from p-hyp", iPhotCount);	
            // verify that rPhoton > rMinKaon(@ phiEstimated = phiPhoton)
            Printf("\npopulate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinProton, getMinCkovProton());");
        isMinKaonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinKaon, getMinCkovKaon(), "Kaon");


            // this means rProtonMin < rPhoton < rMaxPion
            if(isMinKaonOk) {
            
            
                // denne var feil!!! her var det satt arrRMinPion og gtckovMin
                isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");

                if(isMaxPionOk) {

                    isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion");
                    
                //isMinPionOk = ...
                    // this means rPionMin < rPhoton
                    if(isMinPionOk) {
                        isPhotonPionCand = true;
                        pionCandidates.push_back(std::array<double,2>{x,y});	
                        Printf("Photon%d is Pion Candiate", iPhotCount); 
                        // we have pion-candiate
                    }	

                    isMaxKaonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxKaon, getMaxCkovKaon(), "Kaon");
                    // isMaxKaonOk = populate2Ptr->
                    if(isMaxKaonOk) {
                        isPhotonKaonCand = true;
                        kaonCandidates.push_back(std::array<double,2>{x,y});
                        Printf("Photon%d is Kaon Candiate", iPhotCount); 
                        // we have kaon-candiate
                    }							
                }	
                else {
                    Printf("Photon%d could be Pion/Kaon (from p-hyp) but was out of radius-range",iPhotCount);
                }	
            } // end if isMinKaonOk
            else {
                Printf("Photon%d could be Pion/Kaon (from p-hyp) but was out of radius-range",iPhotCount);
            }
        } // end else if (!getProtonStatus() and getKaonStatus() and getPionStatus())



        // denne ser litt rar ut?

        // if neither Proton or Kaon is possible because momentum-threshold is not exceeded:
        else if(!getProtonStatus() and !getKaonStatus() and getPionStatus()) {
            Printf("Photon%d can be Pion from p-hyp", iPhotCount);	

            Printf("populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhotob, arrMaxPion, getMaxCkovPion());");

            isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");
            
            
            if(isMaxPionOk) {

                Printf("populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMinPion, getMaxCkovPion());",iPhotCount);
                
                bool isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion");
                
                if(isMinPionOk) {
                    // we have pion candidate
                    Printf("Photon%d is Pion Candidate", iPhotCount);	
                    isPhotonPionCand = true; 
                    pionCandidates.push_back(std::array<double,2>{x,y});
                } else {
                    Printf("Photon%d could only have been Pion, but didnt fall within radius-range",iPhotCount);
                }
            } 
            else {
                Printf("Photon%d could only have been Pion, but didnt fall within radius-range",iPhotCount);
            }			
            
        } // end else if(!getProtonStatus() and !getKaonStatus() and getPionStatus())

        // this should not really be "possible", but can be if we have another candidate that is not pion/kaon/proton
        else {
            Printf("Photon%d was not Pion Kaon or Proton from p-hyp",iPhotCount);
            // we got particle that should not be able to be pion, kaon or proton
        }


        
        Printf("\n===============================================================");	
        Printf("	Photon%d  Candidate: Pion = %d Kaon = %d Proton %d", iPhotCount, isPhotonPionCand, isPhotonKaonCand, isPhotonProtonCand);
        Printf("===============================================================\n");	



        /* ef: kan vi gjoere noe ala dette for clusters?
        // correctly identified cherenkov photons
        if((isPhotonPionCand or isPhotonKaonCand or isPhotonProtonCand) and etaC != 0) {
        numFoundActualCkov ++;
        }
        
        // correctly identified cherenkov photons
        if(!(isPhotonPionCand or isPhotonKaonCand or isPhotonProtonCand) and etaC != 0) {
        Printf("Ckov %.2f", etaC);
        
        print = true;
            //throw std::invalid_argument("Photon Ckov not found!??");
        } if(etaC != 0) {
            numActualCkov ++;
        }
        // number of found bg photons that are labelled as within proper range:
        if((isPhotonPionCand or isPhotonKaonCand or isPhotonProtonCand) and etaC == 0) {
        numBackgroundLabeledCkov ++;	
        }*/


        // struct with x, y, phiLocal, phi, R + bool isCandidate? 
        // or just x, y, /* phiLocal,*/ phi, R; where x, y... = 0 if !isCandidate

        // better to not store rPhot and phiPhot? this can easily be obtained in ML photon
        /*pionCands.emplace_back(Candidate{x, y, rPhoton, phiPhoton, isPhotonPionCand});	
        kaonCands.emplace_back(Candidate{x, y, rPhoton, phiPhoton, isPhotonKaonCand});	
        protonCands.emplace_back(Candidate{x, y, rPhoton, phiPhoton, isPhotonProtonCand});
        */ 


        int cStatus = 4*static_cast<int>(isPhotonPionCand) + 2*static_cast<int>(isPhotonKaonCand) + 1*static_cast<int>(isPhotonProtonCand);

        //CLusterCandidate :addCandidateStatus(int iTrack, int hadronCandidateBit)
        photons.addCandidateStatus(trackIndex, cStatus);

        /// lagre denne istedet :
        /*
        if( x > 0 && x < 156 && y < 144 && y > 0) {
            candidatesCombined.emplace_back(ParticleUtils::Candidate2{x, y, cStatus});
        } */


        /* ef : TODO check this to verify
        if(true) {
            // TODO: plot as same in maxPion..Range
            if( x > 0 && x < 156 && y < 144 && y > 0) {
                // this means it falls within range
                if(cStatus != 0) {

                    // this means  candidate is a ckov-photon
                    if(etaC != 0) {
                        //(mArrAndMap->ckovCandMapRange)->Fill(x,y);
                        mArrAndMap.fillCkovCandMapRange(x,y); 
                    }	else { // falls within range, but is bg
                    //(mArrAndMap->bgCandMapRange)->Fill(x,y);
                        mArrAndMap.fillbgCandMapRange(x,y);
                    }
                }
                // falls out of range
                else {
                    // this means  candidate is a ckov-photon, but out of range
                    if(etaC != 0) {
                        mArrAndMap.fillckovCandMapOutRange(x,y);
                        //(mArrAndMap->ckovCandMapOutRange)->Fill(x,y);
                    }	else { // falls out of range, and is bg
                        mArrAndMap.fillbgCandMapOutRange(x,y);
                        //(mArrAndMap->bgCandMapOutRange)->Fill(x,y);			
                    }
                }
            }
        } */ 

            // or store as 1 vector, where candidate status 8 => 2^3 (000, 001, 010, 011, 100, 110, 101, 111)
            
            // phiL, phi, R of maxPionVec vectorh



        /*
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        ***************************************************************************
        */ 


            // *<init ClosedForm segm>//		
        // "Closed form approch" init    

        // check if inside pionMax and outside protonMin



            /*
            Printf("\n\n");
        TVector2 rPosPion;
        bool pionBelow = populate.checkRangeBelow(posPhoton, getMaxCkovPion(), rPosPion);

            TVector2 rPosPionMin;
        populate.checkRangeBelow(posPhoton, getMinCkovPion(), rPosPionMin);

        mapPhotons->Fill(posPhoton.X(), posPhoton.Y());
        mapPionMax->Fill(rPosPion.X(), rPosPion.Y());
        mapPionMin->Fill(rPosPionMin.X(), rPosPionMin.Y());

            /	*
        mapPionMaxRev->Fill(rPosPion.X(), rPosPion.Y());
        mapPionMinRev->Fill(rPosPionMin.X(), rPosPionMin.Y());
        * / 


        TVector2 rPosKaonMin, rPosKaonMax;
            populate.checkRangeBelow(posPhoton, getMinCkovKaon(), rPosKaonMin);
            populate.checkRangeBelow(posPhoton, getMaxCkovKaon(), rPosKaonMax);
        mapKaonMax->Fill(rPosKaonMax.X(), rPosKaonMax.Y());
        mapKaonMin->Fill(rPosKaonMin.X(), rPosKaonMin.Y());


        TVector2 rPosProtonMin, rPosProtonMax;
            populate.checkRangeBelow(posPhoton, getMinCkovProton(), rPosProtonMin);
            populate.checkRangeBelow(posPhoton, getMaxCkovProton(), rPosProtonMax);
        mapProtonMax->Fill(rPosProtonMax.X(), rPosProtonMax.Y());
        mapProtonMin->Fill(rPosProtonMin.X(), rPosProtonMin.Y());


        TVector2 rPosProton, temp;
        bool protonBelow = populate.checkRangeBelow(posPhoton, getMinCkovProton(), rPosProton);




        // if not inside ckovMaxPion, continue loop
        if(!pionBelow){
            Printf("!pionBelow");
            continue;
        } if(protonBelow == false && getProtonStatus() == true){
            Printf("protonBelow = 0 && getProtonStatus() = 1");
            continue;
        }
        // global bounds ok, can check candidates 
        else {
            // check Pion
            if(getPionStatus()){ // shouldt really be possible in this case but.. 

                Printf("\n\n ------------ \n\n Region: minProton < ckov < maxProton ok");
                TVector2 rPosPionB;
                bool pionBelow = populate.checkRangeBelow(posPhoton, getMaxCkovPion(), rPosPionB);

                Printf("getPionStatus Ok --> checkRangeAbove getMinCkovPion %.2f ", getMinCkovPion());
                if(populate.checkRangeAbove(posPhoton, getMinCkovPion(), rPosPionB)) {
                mapPion->Fill(xG, yG);
                // add candidate to pions 
                pionCandidates.push_back(std::make_pair(xG, yG));
                Printf("pionCand found ");
                filledBins.push_back(std::make_pair(xG, yG));
                hSignalMap->Fill(xG, yG);
                // range is ok, can check other candidates
                // check if ckov > ckovMaxPion
                }
            }
            if(getKaonStatus()){ // shouldt really be possible in this case but..   		    


            Printf("\n getKaonStatus Ok --> checkRange2 getMinCkovKaon %.2f, getMaxCkovKaon %.2f ", getMaxCkovKaon());

            if (populate.checkRange2(posPhoton, getMaxCkovKaon(), getMinCkovKaon(), temp))
            {	
                Printf("kaonCand found ");			    
                // kaon range ok:
                kaonCandidates.push_back(std::make_pair(xG, yG));
                mapKaon->Fill(xG, yG);
            }
            } 
            if (getProtonStatus()) {
        Printf("\n getProtonStatus Ok --> checkRangeAbove getMaxCkovProton %.2f ", getMinCkovPion());
            if(populate.checkRangeBelow(posPhoton, getMaxCkovProton(), temp)) {
                Printf("protonCand found ");			    
                // proton range ok:
                mapProton->Fill(xG, yG);
                protonCandidates.push_back(std::make_pair(xG, yG));
            }
            } 
            }  // end else / if pionBelow 
            Printf("\n\n");   
                    //<end ClosedForm segm>//		
                    */

        }// end else ifTrue
    }   // end for ckovPhotons


	// iterate through photons to check : 


    /*
    Printf("	Identified the following  : \n");
    Printf("	numBackgroundPhotons %d", numBackgroundPhotons);
    Printf("	numFoundActualCkov %d", numFoundActualCkov);
    Printf("	numActualCkov %d", numActualCkov);
    Printf("	numBackgroundLabeledCkov %d", numBackgroundLabeledCkov);        
    */ 
        
        /*
    if(numFoundActualCkov != numActualCkov) {
                throw std::invalid_argument("wtf value does numActualCkov have?");
    }*/ 

	Printf("=========================================================================");
	Printf("CkovHyps, possible candidates from Momentum | Pion = %d, Kaon = %d, Proton = %d", getPionStatus(), getKaonStatus(), getProtonStatus());


    /*int cntTemp = 0;
  for(const auto& c : candidatesCombined){
		//Printf("Phot%d : x = %.2f, y = %.2f || statusCand = %d", cntTemp++, c.x, c.y, c.candStatus); 
	}*/
	Printf("=========================================================================");

segment
    Printf("number of candidates : proton %zu, kaon %zu, pion %zu", protonCandidates.size(), kaonCandidates.size(), pionCandidates.size());
 
  
    //for(const auto& pair: filledBins)
    //	Printf("CkovTools segment candidates: x%f y%f", pair.first, pair.second);    


    // hNoiseMap->SetMarkerColor(kRed);
    Printf("CkovTools segment filledBins Size %zu", filledBins.size());



    // get impact points of track RAD and PC
    //const auto trkPC = populate.getPcImp();
    /*const auto trkRad = populate.getTrackPos();
    TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
    TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
    trkRadMap->Fill(trkRad.X(), trkRad.Y());
    trkPCMap->Fill(trkPC.X(), trkPC.Y());*/





    /*
    mapProton->SetMarkerColor(kGreen + 4);
    localRefMIP2->Draw("same");  

    TCanvas *segm = new TCanvas("semg","semg",800,800);  
    segm->Divide(2,2);

    segm->cd(1);
    mapPhotons->SetMarkerStyle(2); 

    mapPhotons->Draw();
    localRefMIP2->Draw("same");  
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 


    segm->cd(2);
    mapPion->SetTitle(Form("mapPion : min%d max%d photons%d", mapPionMin->GetEntries(), mapPionMax->GetEntries(), mapPhotons->GetEntries()));
    mapPion->Draw();
    mapPhotons->Draw("same");
    localRefMIP2->Draw("same");
    mapPionMin->Draw("same");
    mapPionMax->Draw("same");
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

    segm->cd(3);

    mapKaon->SetTitle(Form("mapKaon: min%d max%d photons%d", mapKaonMin->GetEntries(), mapKaonMax->GetEntries(), mapPhotons->GetEntries()));
    mapKaon->Draw();
    mapPhotons->Draw("same");
    mapKaonMin->Draw("same");
    mapKaonMax->Draw("same");
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

    segm->cd(4);
    mapProton->SetTitle(Form("mapProton : min%d max%d photons%d", mapProtonMin->GetEntries(), mapProtonMax->GetEntries(), mapPhotons->GetEntries()));
    mapProton->Draw();

    mapPhotons->Draw("same");
    mapProtonMin->Draw("same");
    mapProtonMax->Draw("same");
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

    gPad->Update();
    segm->Show();



    TCanvas *thLocal = new TCanvas("thLocal","thLocal",800,800);  
    thLocal->cd();
    //localBox->Draw();
    //localRef->SetMarkerColor(kBlue);


    localRefMIP->SetMarkerStyle(3);
    localRef->SetMarkerStyle(2);
    localRefUnrot->SetMarkerStyle(2);
    localRefUnrot->SetMarkerColor(kBlue);
    localRefBG->SetMarkerColor(kRed);
    localRefBG->SetMarkerStyle(2);
    localRefMIPUnrot->SetMarkerColor(kBlue);
    localRefMIPUnrot->SetMarkerStyle(3);
    localRefMIPUnrot->Draw("same");
    
    localRef->Draw();
    localRefBG->Draw("same");
    localRefUnrot->Draw("same");

    Printf("ckovtools segment : localPion->Draw()");
    //localPion->Draw("same");    
    localRefMIP->Draw("same");// *
    tlineUpLocal->Draw();
    tlineDownLocal->Draw();
    tlineUpLocalR->Draw();
    tlineDownLocalR->Draw();//  /
    Printf("ckovtools segment : localPion->Draw()");
    gPad->Update();
    // * / 
    */


	if(print) {

		
 		mArrAndMap.setEventCount(eventCnt);
		
		
		Populate2* pPtr = new Populate2(trkPos, trkDir, nF, L);
		mArrAndMap.setPopulatePtr(pPtr);
		
		if(false) { Printf("populate2Ptr was nullptr!");}
	  else {	
		mArrAndMap.setMinArrays(arrMinPionPos, arrMinKaonPos, arrMinProtonPos);
		mArrAndMap.setMaxArrays(arrMaxPionPos, arrMaxKaonPos, arrMaxProtonPos);		
	
	  
	  
		auto d = (rW- L + tGap + qW);


		auto tanP = TMath::Tan(thetaP);
		
		auto cosP = TMath::Cos(phiP);
		auto sinP = TMath::Sin(phiP);
	
		Printf("delta %.2f | x %.2f y %.2f ", d*tanP, d*tanP*cosP, d*tanP*sinP);
		//mArrAndMap->drawTotalMapAndMaxRegions();
		mArrAndMap.drawTotalMap();
		//mArrAndMap->drawMaxRegions();		  
  
		Printf("Event Number%d, momentum %.2f | mass %.2f | thetaP %.2f | phiP %.2f | L %.2f", eventCnt,momentum, mass, thetaP, phiP, L);
					
		// const double& ckovThe, const double& ckovPhi, const double & L
		const auto l1 = populatePtr->tracePhot(getMaxCkovPion(), 0, L);
		const auto l2 = populatePtr->tracePhot(getMinCkovPion(), 3.14159285, L);

		Printf("phiP %.5f: angle rad --> pc = %.5f | Acos %.5f Asin %.5f"  , phiP, (trkPos-trkPC).Phi(), TMath::ACos(d*tanP*cosP/(trkPos-trkPC).Mod()), TMath::ASin(d*tanP*sinP/(trkPos-trkPC).Mod()));
		
		Printf("PHI : l1 %.5f | l2 %.5f ", (l1-trkPos).Phi(),(l2-trkPos).Phi());
		Printf("RADIUS : l1 %.5f | l2 %.5f ", (l1-trkPC).Mod(),(l2-trkPC).Mod());

		Printf("CkovHyps %.2f %.2f %.2f", ckovHyps[0], ckovHyps[1], ckovHyps[2]);
		//Printf("CkovHyps %.2f %.2f %.2f", ckovHyps[0], ckovHyps[1], ckovHyps[2]);

     
		throw std::invalid_argument("print invoked"); 
    }
  // drawTotalMap / drawMaxRegions
  } else {
    // NB! lagre denne istedet : 

    /*
    candCombined = candidatesCombined; // x, y, cStatus

    protonCands = protonCandidates;//.emplace_back()
    kaonCands = kaonCandidates;//.emplace_back()
    pionCands = pionCandidates;//.emplace_back()
   */ 


    // returner disse ogsaa: kun for debugging :

    //arrayInfo = {numBackgroundPhotons, numFoundActualCkov, numActualCkov, numBackgroundLabeledCkov}; 

   	//delete mArrAndMap;	
   }
	return filledBins;
} // end segment


//r0G = rGL + rLP' ; rLP' = Rz(phi)[rLP]
void phiRing2Local(double &xL, double &yL, const double& xMipL, const double& yMipL)
{	  
  TRotation mThetaRot;
  mThetaRot.RotateZ(phiP);
  
  TVector3 mip(xMipL, yMipL, 0);
  TVector3 phiRingPos(xL, yL, 0);
  TVector3 op = mip + mThetaRot*phiRingPos;
  xL = op.Px();
  yL = op.Py();
}
  
//r0G = rGL + rLP' ; rLP' = Rz(phi)[rLP]
// rLP = Rz(-phi) * [r0G-rGL]
void local2PhiRing(double &xL, double &yL, const double& xMipL, const double& yMipL)
{	  
  TRotation mThetaRot;
  mThetaRot.RotateZ(-phiP);
  
  TVector3 pos(xL-xMipL, yL-yMipL, 0);
  TVector3 op = mThetaRot*pos;
  xL = op.Px();
  yL = op.Py();
}

void local2GlobalRef(double& xL, double& yL)
{	  
  xL = xL  + xRad;
  yL = yL  + yRad;	  
}


std::pair<double, double> local2Global(double xL, double yL)
{
  
  const auto x = xL  + xRad;
  const auto y = yL  + yRad;	  
  return {x, y};
}

std::pair<double, double> global2Local(double xG, double yG)
{

  //mTheta.RotateY
  
  const auto x =  xG - xRad;
  const auto y =  yG - yRad;
  return std::make_pair(x, y);
}




// NB! not uses L by simulation here!
// get radius from MIP to a specific phiL, etaC pair
double getR(double etaC, double phiL)
{

	const auto cosPhiL = TMath::Cos(phiL); // 
	const auto sinPhiL = TMath::Cos(phiL); // --||-- 
	
	const auto cosEtaC = TMath::Cos(etaC);
	const auto sinEtaC = TMath::Sin(etaC);

	const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));
	const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));
	const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC * (rwlDeltaR + qwDeltaR));


 // ef :error was on this line :
 // 		const auto denum = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto tZ = num/denum;

	const auto Lz = (rW-L) + qW + tZ;

	const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
  //Printf("getR_Lmax : R %f |  wlGap %f qwDeltaR %f tGapDeltaR %f", R, rwlDeltaR,qwDeltaR,tGapDeltaR);
	return R;
}

	// get radius from MIP to a specific phiL, etaC pair
double getR_Lmax(double etaC, double phiL)
{ 


	const auto cosPhiL = TMath::Cos(phiL); // 
	const auto sinPhiL = TMath::Cos(phiL); // --||-- 
	
	const auto cosEtaC = TMath::Cos(etaC);
	const auto sinEtaC = TMath::Sin(etaC);

	const auto rwlDeltaR = (rW - lMax)/(TMath::Sqrt(1-sinEtaC*sinEtaC));
	const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));
	const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC * (rwlDeltaR + qwDeltaR));

	const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));
	const auto tZ = num/denum;
	const auto Lz = (rW-lMax) + qW + tZ;

	const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
  //Printf("getR_Lmax : R %f |  wlGap %f qwDeltaR %f tGapDeltaR %f", R, rwlDeltaR,qwDeltaR,tGapDeltaR);
	return R;
}






// placeholder...
void setArrayMax(double etaTRS, vecArray4& inPutVectorAngle, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();
  const float lMin = 0.;
  const auto trkPC2 = populatePtr->getPcImp();
  const auto trkRad2 = populatePtr->getTrackPos();


  //Printf("\n\n setArrayMax() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = Double_t(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		
		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		//Printf("\n setArrayMax() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS

		// make this fcn in populate instead?		
		populatePtr->trs2Lors(dirTrs, thetaR, phiR);
		
		//Printf("setArrayMax() called trs2Lors on dirTRS : return { thetaR = %.2f, phiR = %.2f}", thetaR, phiR);
		

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		//Printf("setArrayMax() dirLORS {x %.2f y %.2f z %.2f}", dirLORS.X(), dirLORS.Y(), dirLORS.Z());


		// temp
		// ckovThe, const double& ckovPhi, const double & L
		// this should return the same as max
		//const auto t = populatePtr->tracePhot(etaTRS, phiL, lMin);

		// temp
		
		//Printf("setArrayMax() called  populatePtr->traceForward(dirLORS (thetaR %.2f, phiR %.2f) lMin =  %.2f", thetaR, phiR, lMin);
		const auto& max = populatePtr->traceForward(dirLORS, lMin);
		
		const auto r = (max - trkPC2).Mod();
		auto phiPC = (max - trkPC2).Phi();
		//Printf("setArrayMax() traceForward returned TVector2 {x %.2f y %.2f} == > R  = %.2f", max.X(), max.Y(), r);
		//Printf("setArrayMax() tracePhot returned TVector2 {x %.2f y %.2f} == > R  = %.2f", t.X(), t.Y(), (t-trkPC2).Mod());
		
		// add protection if traceForward returns 0 or -999?


 		//Printf("setArrayMax() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC.X(), trkPC.Y());
 		

		// check at phiR  er det samme som (t-rad).phi og (max-rad).phi?
		



		//Printf("setArrayMax2() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC2.X(), trkPC2.Y());


		// phiR in [-pi, pi]? set to 0..2pi?
		// inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
		if(phiR < 0) {
			phiR = TMath::TwoPi() + phiR;

    }
     
		if(phiPC < 0) {
			phiPC = TMath::TwoPi() + phiPC;

    }   

    
		/*if( (max-trkRad2).Phi() != (t-trkRad2).Phi()) {
			throw std::invalid_argument("(max-trkRad2).Phi() != (t-trkRad2).Phi())");
		}

		if(TMath::Abs(phiR - (t-trkRad2).Phi()) > 0.0001 && t.X() != -999 && t.Y() != -999)   
		{
		
		  Printf("phiR at %.6f | (t-rad) %.6f | (max-rad) %.6f ", phiR, (t-trkRad2).Phi(), (max-trkRad2).Phi());
			throw std::invalid_argument("phiR != (t-trkRad2).Phi())");
		}
		// protections if r > value?
		//if(r > )*/


		// TODO: if it goes out of map, find intersection with chamber-edges??
		// really jsut check wether TraceForward returned x, y = -999.
 		// to set points out of the map here is fine	


		if((max.Y() == -999) or (max.X() == -999)) {
    	//inPutVector[i] = {0,0,0, 0, 0};
    	//inPutVector.emplace_back(std::array<double, 3>{0, 0, 0});
			// placeholder, find better solution? 
			if(max.Y() == -999) {
				//Printf("setArrayMax() max.Y() %.1f == -999", max.Y());
			}
			if(max.X() == -999) {
				//Printf("setArrayMax() max.X() %.1f == -999", max.X());
			}
    } else {
    	
    	//inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	
    	inPutVectorPos.emplace_back(std::array<double, 2>{max.X(), max.Y()}); 
    	inPutVectorAngle.emplace_back(std::array<double, 4>{phiL, phiR, phiPC, r});   	
    	//inPutVectorPos[i] = {max.X(), max.Y()};
    	//inPutVectorAngle[i] = {phiL, phiR, r};
    	
    	
    	//inPutVectorAngle.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	//Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);
    }
    	//inPutVector[i] = {phiL, phiR, r};
    // Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);

		/*if(maxPion.X() > 0 && maxPion.X() < 156.0 && maxPion.Y() > 0 && maxPion.Y() < 144) {
			// hMaxPion->Fill(maxPion.X(), maxPion.Y());
			// maxPionVec[i] = std::make_pair(maxPion.X(), maxPion.Y());
			Printf("maxPion loop i = %d, maxSize = %zu", i, maxPionVec.size()); 
		}*/
	}/*
  	Printf("\n");
	for(const auto& ip : inPutVectorAngle) {
		const auto& phiL_ = ip[0];
		const auto& phiR_ = ip[1]; 
		const auto& r_ = ip[2];  
		//Printf("setArrayMax() --> checking inputVector | : phiL %.2f, phiR %.2f, r %.2f", phiL_, phiR_, r_);
		if(r_ == 0 ) {throw std::invalid_argument("r====??????;");}
  } */
}




// placeholder...
void setArrayMin(double etaTRS, vecArray4& inPutVectorAngle, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();

  // higher L == > lower radius
  const auto lMax = 1.5;
  const auto trkPC2 = populatePtr->getPcImp();

 // Printf("\n\n setArrayMin() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = Double_t(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;

		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		//Printf("\n setArrayMin() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS

		// make this fcn in populate instead?		
		populatePtr->trs2Lors(dirTrs, thetaR, phiR);
		
	  //Printf("setArrayMin() called trs2Lors on dirTRS : return { thetaR = %.2f, phiR = %.2f}", thetaR, phiR);
		

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		//Printf("setArrayMin() dirLORS {x %.2f y %.2f z %.2f}", dirLORS.X(), dirLORS.Y(), dirLORS.Z());


		// temp
		// ckovThe, const double& ckovPhi, const double & L
		// this should return the same as max
		//const auto t = populatePtr->tracePhot(etaTRS, phiL, lMax);
		//Printf("setArrayMin() tracePhot returned TVector2 {x %.2f y %.2f}", t.X(), t.Y());
		// temp
		
		// Printf("setArrayMin() called  populatePtr->traceForward(dirLORS (thetaR %.2f, phiR %.2f) lMax =  %.2f", thetaR, phiR, lMax);
		const auto& max = populatePtr->traceForward(dirLORS, lMax);
		
		// Printf("setArrayMin() traceForward returned TVector2 {x %.2f y %.2f}", max.X(), max.Y());
		
		// add protection if traceForward returns 0 or -999?
		const auto r = (max - trkPC2).Mod();

 		// Printf("setArrayMin() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC.X(), trkPC.Y());

		// Printf("setArrayMin2() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC2.X(), trkPC2.Y());


		// phiR in [-pi, pi]? set to 0..2pi?
		// inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
		if(phiR < 0) {
			phiR = TMath::TwoPi() + phiR;
    }
		auto phiPC = (max - trkPC2).Phi();
		if(phiPC < 0) {
			phiPC = TMath::TwoPi() + phiPC;

    }   


		// protections if r > value?
		//if(r > )


		// TODO: if it goes out of map, find intersection with chamber-edges??
		// really jsut check wether TraceForward returned x, y = -999.
 		// to set points out of the map here is fine	


		if((max.Y() == -999) or (max.X() == -999)) {

			// set element to false, make better solution later..
    	//inPutVector[i] = {0, 0, 0};
			// placeholder, find better solution? s
			if(max.Y() == -999) {
				// sPrintf("setArrayMin() max.Y() %.1f == -999", max.Y());
			}
			if(max.X() == -999) {
				// Printf("setArrayMin() max.X() %.1f == -999", max.X());
			}
    } else {
    	inPutVectorPos.emplace_back(std::array<double, 2>{max.X(), max.Y()}); 
    	inPutVectorAngle.emplace_back(std::array<double, 4>{phiL, phiR, phiPC, r});   	
    	//inPutVectorPos[i] = {max.X(), max.Y()};
    	//inPutVectorAngle[i] = {phiL, phiR, r};
    	//Printf("setArrayMin() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);
    }
    	//inPutVector[i] = {phiL, phiR, r};
    // Printf("setArrayMin() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);

		/*if(maxPion.X() > 0 && maxPion.X() < 156.0 && maxPion.Y() > 0 && maxPion.Y() < 144) {
			// hMaxPion->Fill(maxPion.X(), maxPion.Y());
			// maxPionVec[i] = std::make_pair(maxPion.X(), maxPion.Y());
			Printf("maxPion loop i = %d, maxSize = %zu", i, maxPionVec.size()); 
		}*/
	}/*
  	Printf("\n");
	for(const auto& ip : inPutVectorAngle) {
		const auto& phiL_ = ip[0];
		const auto& phiR_ = ip[1]; 
		const auto& r_ = ip[2];  
		 Printf("setArrayMin() --> checking inputVector | : phiL %.2f, phiR %.2f, r %.2f", phiL_, phiR_, r_);
  } */ 
}


// map i.e., hMaxPion | vecArray i.e. arrMaxPion
 

// placeholder...
void setArrayMin2(Populate* populate, double etaTRS, vecArray3& inPutVector)
{
  const size_t kN = inPutVector.size();
  const auto lMax = 1.5;
		
	for(int i = 0; i < kN; i++){

		const auto phiL = Double_t(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);
		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS
		populatePtr->trs2Lors(dirTrs, thetaR, phiR);
		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		const auto& max = populatePtr->traceForward(dirLORS, lMax);

		const auto r = (max - trkPC).Mod();
		inPutVector[i] = {phiL, phiR, r};


	}
}


void populateRegions(std::vector<std::pair<double, double>>& vecArr, TH2F* map, const double& eta, const double& l) {	 
	 const int kN = vecArr.size();
	 for(int i = 0; i < vecArr.size(); i++){
		  const auto& value = populatePtr->tracePhot(eta, Double_t(TMath::TwoPi()*(i+1)/kN), l);
		  if(/*value.X() > 0 && value.X() < 156.0 && value.Y() > 0 && value.Y() < 144*/true) {
		    map->Fill(value.X(), value.Y());
		    vecArr[i] = std::make_pair(value.X(), value.Y());
		  }
	 }
 		
 }


}; // end class CkovTools






