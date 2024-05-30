// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
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
/*  virtual void process(gsl::span<o2::tpc::Digit const> const& digits,
   o2::dataformats::ConstMCLabelContainerView const& mcDigitTruth) = 0;
 */
//,
// void Clusterer::Dig2Clu(gsl::span<const o2::hmpid::Digit> digs, std::vector<o2::hmpid::Cluster>& clus, float* pUserCut, bool isUnfold, o2::dataformats::ConstMCLabelContainerView const& mcDigitTruth, )

void Clusterer::Dig2Clu(gsl::span<const o2::hmpid::Digit> digs, std::vector<o2::hmpid::Cluster>& clus, float* pUserCut, MCLabelContainer const* digitMCTruth, bool isUnfold)
{

  // Finds all clusters for a given digits list provided not empty. Currently digits list is a list of all digits for a single chamber.
  // Puts all found clusters in separate lists, one per clusters.
  // Arguments: pDigAll     - list of digits for all chambers
  //            pCluAll     - list of clusters for all chambers
  //            isTryUnfold - flag to choose between CoG and Mathieson fitting
  //  Returns: none

  const int numCluStart = clus.size();


  /*if (digitMCTruth == nullptr) {
    LOGP(info, "digitMCTruth was nullptr!");
  }*/

  struct Pad {
    int x, y, m;
  };

  TMatrixF padMap(Param::kMinPx, Param::kMaxPcx, Param::kMinPy, Param::kMaxPcy); // pads map for single chamber 0..159 x 0..143
  int pUsedDig = -1;
  int padChX = 0, padChY = 0, module = 0;
  std::vector<Pad> vPad;
  std::vector<const Digit*> digVec;
  for (int iCh = Param::kMinCh; iCh <= Param::kMaxCh; iCh++)

  { // chambers loop

    padMap = (Float_t)-1; // reset map to -1 (means no digit for this pad)

    for (size_t iDig = 0; iDig < digs.size(); iDig++) {
      o2::hmpid::Digit::pad2Absolute(digs[iDig].getPadID(), &module, &padChX, &padChY);
      vPad.push_back({padChX, padChY, module});

      if (module == iCh) {
        padMap(padChX, padChY) = iDig; // fill the map for the given chamber, (padx,pady) cell takes digit index
      }

    } // digits loop for current chamber


    for (size_t iDig = 0; iDig < digs.size(); iDig++) { // digits loop for current chamber

      if (vPad.at(iDig).m != iCh || (pUsedDig = UseDig(vPad.at(iDig).x, vPad.at(iDig).y, padMap)) == -1) { // this digit is from other module or already taken in FormClu(), go after next digit
        continue;
      }

      digVec.clear();
      Cluster clu;
      clu.setDigits(&digVec);
      clu.setCh(iCh);

      std::vector<int> digitIndicesRawCluster;
      // ef : added, for MC > digit indices of raw cluster (before solving)      
      
      // if using MC 
      if (digitMCTruth != nullptr) {
        // pass indices by reference, keep track of indices of digits used in cluster
        FormCluMC(clu, pUsedDig, digs, padMap, digitIndicesRawCluster); // form cluster starting from this digit by recursion
      } 
      // not using MC
      else {
        FormClu(clu, pUsedDig, digs, padMap); // form cluster starting from this digit by recursion
      }

      // recursively add all adjacent digits to the cluster
      int formedClusters = -1;
      const int cluSizeIn = clus.size();

      if (digitMCTruth == nullptr) {
        formedClusters = clu.solve(&clus, pUserCut, isUnfold); // solve this cluster and add all unfolded clusters to provided list
      }

      // for setting cluster MC-labels
      if (digitMCTruth != nullptr) {

        std::map<int, std::vector<int>> resolvedDigIndicesMap;

        // solving convoluted clusters, and getting local indices
        //  each entry in the resolvedDigIndicesMap is for each resolved cluster a vector
        // containing the indices of the their resolved digits (the digits forming the resolved cluster)

        formedClusters = clu.solveMC(&clus, pUserCut, isUnfold, resolvedDigIndicesMap); // solve this cluster and add all unfolded clusters to provided list

        const int startIndex = cluSizeIn;

        // if (formedClusters > 1 && formedClusters < 6) {

        // we loop over the number of formed clusters(formedClusters)
        // and set the MC labels per resoolved clusteers
        if (formedClusters > 1 && formedClusters < kMaxLocMax) {

          // use resolvedDigIndicesMap
          // for(int i = 0; i < resolvedClusters; i++) {
          //  cosnt auto& cluster = clus[i+cluSizeIn];
          //  vector<int> resolvedDigIndices = resolvedDigIndicesMap[i];
          //
          // map from "local" to global



          // loop over all the resolved clusters
          // resolvedDigIndices is the (int, vector) combination
          //  resolvedDigIndices.first is the index of the resolved cluster
          // resolvedDigIndices.second is the vector of the indices of their resolved digits
          for (const auto& resolvedDigIndices : resolvedDigIndicesMap) {
            // 0.1.2....

            // the integer in the map;
            // gives the index of the current resolved index out of the formed clusters (formedClusters)
            int indSolvedClu = resolvedDigIndices.first;

            // get cluster
            // current entry is the end of current cluster-vector
            // and adding the index of the resolved cluster out of the numResolvedClusters
            int currentEntry = startIndex + indSolvedClu;


            // ef > todo : check on size
            const auto& cluster = clus[currentEntry];

            std::vector<int> resolvedDigIndsVec = resolvedDigIndices.second;
            // get all the digits, and all the MC-labels from them


            // iterate over the local to map to global and fidn selected
            // store the global indicwes of the resolved digits for the current resolved cluster
            std::vector<int> globalIndResolvedDigits;
            for (const auto& index : resolvedDigIndsVec) {
              globalIndResolvedDigits.push_back(digitIndicesRawCluster[index]);
            }

            if (mClsLabels) {
              addCluLabelsFromDig(cluster, globalIndResolvedDigits, digitMCTruth, mClsLabels, clus.size());

            }
          }


        } // end if (formedClusters > 1 && formedClusters < 6) {

        // in this case, only 1 cluster is formed
        // because kMaxLocMax was exceeded
        // in this case we use the raw cluster
        else if (formedClusters == 1 || formedClusters >= kMaxLocMax) { // 6 is max allowed value kMaxLocMax

          //  const auto& cluster = clus[i+cluSizeIn];
          // the "unresolved" are here also the resolved,
          // since we only had 1 cluster to begin with in raw cluster
          // we can directly use digitIndicesRawCluster as the global
          const auto& cluster = clus[startIndex];
          addCluLabelsFromDig(cluster, digitIndicesRawCluster, digitMCTruth, mClsLabels, clus.size());

        } // end else if (formedClusters == 1 || formedClusters >= 6 ) {

      } // if digitsMcTruth

    } // digits loop for current chamber
    vPad.clear();
  } // chambers loop

  // const int numCluStart = clus.size();


  startIndexDigMC += digs.size();

  return;

} // Dig2Clu()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Clusterer::addCluLabelsFromDig(const Cluster& cluster, const std::vector<int>& globalDigitIndices, MCLabelContainer const* digitMCTruth, MCLabelContainer* mClsLabels, int cluSize)
{

  int lbl = mClsLabels->getIndexedSize(); // this should correspond to the current number of clusters

  // loop over resolvedIndexes array;
  int numDigitsCount = 0;
  for (int i = 0; i < globalDigitIndices.size(); i++) {
    const int& indexOfDigGlobal = globalDigitIndices[i];

    numDigitsCount++;
    // const auto& digOfClu = &digs[iDig];

    if (digitMCTruth == nullptr) {
      continue;
    }

    int digitLabel = startIndexDigMC + indexOfDigGlobal;
    // all MC-hits from the specific digit
    gsl::span<const o2::MCCompLabel> mcArray = digitMCTruth->getLabels(digitLabel);


    for (int j = 0; j < static_cast<int>(mcArray.size()); j++) {

      const auto& currentIndex = digitMCTruth->getMCTruthHeader(digitLabel).index + j;
      auto label = digitMCTruth->getElement(currentIndex);

      // same as digits having multiple hits
      // clsuters have multiple digits
      // we fill MC-Complabel label at index lbl for headArray
      // this is the hit for a digit in the cluster

      mClsLabels->addElement(lbl, label);

    } // end for mcArray
  }   // end for iDig : resolvedDigIndices


} // end addCluLabelsFromDig
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

  } // neighbours loop

} // FormClu()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ef : added function>
// for MC labeling, the propagation of MC labels of the digits forming the cluster is done here
void Clusterer::FormCluMC(Cluster& pClu, int pDig, gsl::span<const o2::hmpid::Digit> digs, TMatrixF& pDigMap, std::vector<int>& digitIndicesRawCluster)
{
  // Forms the initial cluster as a combination of all adjascent digits. Starts from the given digit then calls itself recursevly  for all neighbours.
  // Arguments: pClu - pointer to cluster being formed
  //   Returns: none
  //   void digAdd(const o2::hmpid::Digit* pDig);                                         // add new digit to the cluster
  pClu.digAdd(&digs[pDig]); // take this digit in cluster
  digitIndicesRawCluster.push_back(pDig);
  // ef : add index of digit to cluster
  // pClu.setUnresolvedIndex(pDig);

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
  }                               // up
  for (int i = 0; i < cnt; i++) { // neighbours loop
    pDig = UseDig(cx[i], cy[i], pDigMap);
    if (pDig != -1) {
      FormCluMC(pClu, pDig, digs, pDigMap, digitIndicesRawCluster);
    } // check if this neighbour pad fired and mark it as taken
  }   // neighbours loop
} // FormCluMC()

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
