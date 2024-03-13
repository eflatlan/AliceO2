// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Clusterer.cxx
/// \brief Implementation of the HMPID cluster finder
#include <algorithm>
#include "FairLogger.h" // for LOG
#include "Framework/Logger.h"
#include "HMPIDBase/Param.h"
#include "DataFormatsHMP/Cluster.h"
#include "HMPIDReconstruction/Clusterer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include <TStopwatch.h>

using namespace o2::hmpid;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*Clusterer::Clusterer()
{
  // mDigs = {nullptr};
  // mClus = {nullptr};
}*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*  virtual void process(gsl::span<o2::tpc::Digit const> const& digits, o2::dataformats::ConstMCLabelContainerView const& mcDigitTruth) = 0;
*/
//, 
//void Clusterer::Dig2Clu(gsl::span<const o2::hmpid::Digit> digs, std::vector<o2::hmpid::Cluster>& clus, float* pUserCut, bool isUnfold, o2::dataformats::ConstMCLabelContainerView const& mcDigitTruth, )
void Clusterer::Dig2Clu(gsl::span<const o2::hmpid::Digit> digs, std::vector<o2::hmpid::Cluster>& clus, float* pUserCut, MCLabelContainer const* digitMCTruth, bool isUnfold)
{

  // Finds all clusters for a given digits list provided not empty. Currently digits list is a list of all digits for a single chamber.
  // Puts all found clusters in separate lists, one per clusters.
  // Arguments: pDigAll     - list of digits for all chambers
  //            pCluAll     - list of clusters for all chambers
  //            isTryUnfold - flag to choose between CoG and Mathieson fitting
  //  Returns: none

  LOGP(info, "\n\n called Dig2Clu : ");
  LOGP(info, " clus Size {}", clus.size());
  LOGP(info, " digs Size {}", digs.size());  
  
  if(digitMCTruth == nullptr) {
    LOGP(info, "digitMCTruth was nullptr!");
  }
  
  struct Pad {
    int x, y, m;
  };

  TMatrixF padMap(Param::kMinPx, Param::kMaxPcx, Param::kMinPy, Param::kMaxPcy); // pads map for single chamber 0..159 x 0..143

  int pUsedDig = -1;
  int padChX = 0, padChY = 0, module = 0;
  std::vector<Pad> vPad;
  std::vector<const Digit*> digVec;
  for (int iCh = Param::kMinCh; iCh <= Param::kMaxCh; iCh++) { // chambers loop
    padMap = (Float_t)-1;                                      // reset map to -1 (means no digit for this pad)
    for (size_t iDig = 0; iDig < digs.size(); iDig++) {
      o2::hmpid::Digit::pad2Absolute(digs[iDig].getPadID(), &module, &padChX, &padChY);
      vPad.push_back({padChX, padChY, module});
      if (module == iCh) {
        padMap(padChX, padChY) = iDig; // fill the map for the given chamber, (padx,pady) cell takes digit index
      }
    } // digits loop for current chamber
    LOGP(info, "looping over digits");
    for (size_t iDig = 0; iDig < digs.size(); iDig++) { // digits loop for current chamber
      // o2::hmpid::Digit::pad2Absolute(digs[iDig].getPadID(), &module, &padChX, &padChY);
      if (vPad.at(iDig).m != iCh || (pUsedDig = UseDig(vPad.at(iDig).x, vPad.at(iDig).y, padMap)) == -1) { // this digit is from other module or already taken in FormClu(), go after next digit
        continue;
      }
      digVec.clear();
      Cluster clu;
      clu.setDigits(&digVec);


      // denne skal mulgiens vaere nptr, men jeg skjonner ikke logikken!
      /*if(!clu.dig(0))
    	Printf("dig2Clu Digit i0 nullptr ");
      else 
      	Printf("dig2Clu Digit i0 ok"); */

      clu.setCh(iCh);


      //ef : ad digitMCTruth
      FormClu(clu, pUsedDig, digs, padMap); // form cluster starting from this digit by recursion



      // here we should propagate digit-MC truth label 
      // selected digits are gotten with pDig->getDigits(); which is a subset of Digits from digs; // then we should iterate over those to setMCTruth... 
      /*
      if(digitMCTruth != nullptr) {
        clu.setMCTruth(digitMCTruth); // we do this in Cluster.h
      } */

      
      // filling the MC labels of this cluster; the first will be those of the main digit; then the others
      if (digitMCTruth != nullptr) {
          int lbl = mClsLabels->getIndexedSize(); // this should correspond to the number of digits also;


          LOGP(info, "\n================================\n");

          LOGP(info, "lbl = {} (mClsLabels->getIndexedSize())", lbl);



          const auto& digsOfClu = *(clu.getDigits());

          for(int digIndex = 0; digIndex < digsOfClu.size(); digIndex++) {  
              LOGP(info, "contributing digit = {}", digIndex);

              int digitLabel = digsOfClu[digIndex]->getLabel();
              const int digEventNum = digsOfClu[digIndex]->getEventNumber();

              //printf("digitLabel = %d\n", digitLabel);

              LOGP(info, "digitLabel  = {}", digitLabel);

              gsl::span<const o2::MCCompLabel> mcArray = digitMCTruth->getLabels(digitLabel);
              for (int j = 0; j < static_cast<int>(mcArray.size()); j++) {
                auto label = digitMCTruth->getElement(digitMCTruth->getMCTruthHeader(digitLabel).index + j);
                mClsLabels->addElement(lbl, label);
                LOGP(info, "checking element {} in array of labels", j);

                //printf("checking element %d in the array of labels\n", j);
                //printf("EventID = %d\n", label.getEventID());
                LOGP(info, "EventID from MC-label  = {}",  label.getEventID());
                LOGP(info, "EventID from dig:  = {}", digEventNum);

              }
          }
      }
      

      // ,  MCLabelContainer const* digitMCTruth
      // clu now has a vector<Digit*>* field;
      // iterate over it and find cluster topology

      // denne skal ikke vaere nptr

      /*
      if(clu.dig(0)==nullptr)
    	Printf("dig2Clu Digit i0 nullptr ");
      else 
      	Printf("dig2Clu Digit i0 ok");
      */
	
			LOGP(info, "set CluCh {}", iCh);

      clu.setCh(iCh);


      clu.solve(&clus, pUserCut, isUnfold); // solve this cluster and add all unfolded clusters to provided list
			LOGP(info, "clu.solve ok");
      /*
      if(clu.dig(0)==nullptr)
    	Printf("dig2Clu Digit i0 nullptr ");
      else 
      	Printf("dig2Clu Digit i0 ok");
      
      if(clus.back().dig(0) == nullptr) {Printf("dig2Clu dig was nullptr!!");}
      */
    }                                       // digits loop for current chamber
    vPad.clear();
  } // chambers loop
  return;
} // Dig2Clu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*void Clusterer::Dig2Clu(gsl::span<const o2::hmpid::Digit> digs, std::vector<o2::hmpid::Cluster>& clus, float* pUserCut, bool isUnfold)
{
  // Finds all clusters for a given digits list provided not empty. Currently digits list is a list of all digits for a single chamber.
  // Puts all found clusters in separate lists, one per clusters.
  // Arguments: pDigAll     - list of digits for all chambers
  //            pCluAll     - list of clusters for all chambers
  //            isTryUnfold - flag to choose between CoG and Mathieson fitting
  //  Returns: none

  TMatrixF padMap(Param::kMinPx, Param::kMaxPcx, Param::kMinPy, Param::kMaxPcy); // pads map for single chamber 0..159 x 0..143

  int pUsedDig = -1;
  int padChX = 0, padChY = 0, module = 0;

  for (int iCh = Param::kMinCh; iCh <= Param::kMaxCh; iCh++) { // chambers loop
    padMap = (Float_t)-1;                                      // reset map to -1 (means no digit for this pad)
    for (size_t iDig = 0; iDig < digs.size(); iDig++) {
      o2::hmpid::Digit::pad2Absolute(digs[iDig].getPadID(), &module, &padChX, &padChY);
      if (module == iCh) {
        padMap(padChX, padChY) = iDig; // fill the map for the given chamber, (padx,pady) cell takes digit index
      }
    } // digits loop for current chamber

    for (size_t iDig = 0; iDig < digs.size(); iDig++) { // digits loop for current chamber
      o2::hmpid::Digit::pad2Absolute(digs[iDig].getPadID(), &module, &padChX, &padChY);
      if (module != iCh || (pUsedDig = UseDig(padChX, padChY, padMap)) == -1) { // this digit is from other module or already taken in FormClu(), go after next digit
        continue;
      }
      Cluster clu;
      clu.setCh(iCh);
      FormClu(clu, pUsedDig, digs, padMap); // form cluster starting from this digit by recursion
      clu.solve(&clus, pUserCut, isUnfold); // solve this cluster and add all unfolded clusters to provided list
    }                                       // digits loop for current chamber
  }                                         // chambers loop
  return;
} // Dig2Clu()
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Clusterer::FormClu(Cluster& pClu, int pDig, gsl::span<const o2::hmpid::Digit> digs, TMatrixF& pDigMap)
{
  // Forms the initial cluster as a combination of all adjascent digits. Starts from the given digit then calls itself recursevly  for all neighbours.
  // Arguments: pClu - pointer to cluster being formed
  //   Returns: none

  //   void digAdd(const o2::hmpid::Digit* pDig);                                         // add new digit to the cluster
  pClu.digAdd(&digs[pDig]); // take this digit in cluster

  int cnt = 0;
  int cx[4];
  int cy[4];
  int padChX = 0;
  int padChY = 0;
  int module = 0;
  o2::hmpid::Digit::pad2Absolute(digs[pDig].getPadID(), &module, &padChX, &padChY);



  if (padChX > Param::kMinPx) {
    cx[cnt] = padChX - 1;
    cy[cnt] = padChY;
    cnt++;
  } // left
  if (padChX < Param::kMaxPcx) {
    cx[cnt] = padChX + 1;
    cy[cnt] = padChY;
    cnt++;
  } // right
  if (padChY > Param::kMinPy) {
    cx[cnt] = padChX;
    cy[cnt] = padChY - 1;
    cnt++;
  } // down
  if (padChY < Param::kMaxPcy) {
    cx[cnt] = padChX;
    cy[cnt] = padChY + 1;
    cnt++;
  } // up

  for (int i = 0; i < cnt; i++) { // neighbours loop
    pDig = UseDig(cx[i], cy[i], pDigMap);
    if (pDig != -1) {
      FormClu(pClu, pDig, digs, pDigMap);
    } // check if this neighbour pad fired and mark it as taken
  }   // neighbours loop
} // FormClu()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Clusterer::UseDig(int padX, int padY, TMatrixF& pPadMap)
{
  // Digit map contains a matrix if digit numbers.
  // Main operation in forming initial cluster is done here. Requested digit pointer is returned and this digit marked as taken.
  // Arguments: padX,padY - pad number
  //            pDigLst   - list of digits for one sector
  //            pDigMap   - map of those digits
  //   Returns: index to digit if not yet used or -1 if used
  int iDig = (int)pPadMap(padX, padY);
  pPadMap(padX, padY) = -1; // take digit number from the map and reset this map cell to -1
  return iDig;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// bool Clusterer::IsDigSurvive(Int_t *pUserCut, Digit *pDig)const
bool Clusterer::IsDigSurvive(Digit* pDig) const
{
  // Check if the current digit survive to a riapllied sigma cut
  // Arguments: pDig pointer to the current digit
  //   Returns: kTRUE if charge > mean+n*sigma
  /*int iCh = pDig->Ch();
  int iDaqSigCut =(int)fDaqSig->At(iCh)->GetUniqueID();
  if(pUserCut[iCh]<=iDaqSigCut) return kTRUE;
  TMatrixF *pM = (TMatrixF*)fDaqSig->At(pDig->Ch());
  float sig = (*pM)(pDig->PadChX(),pDig->PadChY());
  //if(pDig->Q()>pUserCut[iCh]*sig) return kTRUE; to be improved
  */
  if (pDig->getQ() > 4.) {
    return true;
  } else {
    return false;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*void Clusterer::process(std::vector<o2::hmpid::Digit> const& digits, std::vector<o2::hmpid::Cluster>& clusters, MCLabelContainer const* digitMCTruth)
{
  TStopwatch timerProcess;
  timerProcess.Start();

  //  reader.init();
  //  int totNumDigits = 0;
  //
  //  while (reader.getNextStripData(mStripData)) {
  //    LOG(debug) << "HMPIDClusterer got Strip " << mStripData.stripID << " with Ndigits "
  //               << mStripData.digits.size();
  //    totNumDigits += mStripData.digits.size();
  //
  //    processStrip(clusters, digitMCTruth);
  //  }

  //  LOG(debug) << "We had " << totNumDigits << " digits in this event";
  timerProcess.Stop();
  printf("Timing:\n");
  printf("Clusterer::process:        ");
  timerProcess.Print();
}*/
