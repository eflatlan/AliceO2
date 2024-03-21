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

#include "Steer/MCKinematicsReader.h"

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

  LOGP(info, "\n\n ============================== called Dig2Clu\n looping over digits");

  LOGP(info, " clus Size {}", clus.size());

  LOGP(info, " digs Size {}", digs.size());

  if (digitMCTruth == nullptr)

  {

    LOGP(info, "digitMCTruth was nullptr!");
  }

  struct Pad

  {

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

    for (size_t iDig = 0; iDig < digs.size(); iDig++)

    {

      o2::hmpid::Digit::pad2Absolute(digs[iDig].getPadID(), &module, &padChX, &padChY);

      vPad.push_back({padChX, padChY, module});

      if (module == iCh)

      {

        padMap(padChX, padChY) = iDig; // fill the map for the given chamber, (padx,pady) cell takes digit index
      }

    } // digits loop for current chamber

    LOGP(info, "\n\n ============================== \n looping over digits ============================================================\n ==============================\n");

    for (size_t iDig = 0; iDig < digs.size(); iDig++)

    { // digits loop for current chamber

      // o2::hmpid::Digit::pad2Absolute(digs[iDig].getPadID(), &module, &padChX, &padChY);

      // LOGP(info, "\n\n ================================\n New iDig : {} \n================================\n\n", iDig);

      if (vPad.at(iDig).m != iCh || (pUsedDig = UseDig(vPad.at(iDig).x, vPad.at(iDig).y, padMap)) == -1)

      { // this digit is from other module or already taken in FormClu(), go after next digit

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

      std::vector<int> indicesUnresolved;
      // digit indices of raw cluster (before solving)

      if (digitMCTruth != nullptr) {
        // pass indices by reference, keep track of indices of digits used in cluster
        FormCluMC(clu, pUsedDig, digs, padMap, indicesUnresolved); // form cluster starting from this digit by recursion
      } else {
        FormClu(clu, pUsedDig, digs, padMap); // form cluster starting from this digit by recursion
      }

      // recursively add all adjacent digits to the cluster

      //  ef: added setting the indexes of the digits used in cluster

      // here we should propagate digit-MC truth label

      // selected digits are gotten with pDig->getDigits(); which is a subset of Digits from digs;

      // then we should iterate over those to setMCTruth...

      /*







      if(digitMCTruth != nullptr) {



        clu.setMCTruth(digitMCTruth); // we do this in Cluster.h



      } */

      // filling the MC labels of this cluster; the first will be those of the main digit; then the others

      // will be nullptr if useMc in Digits2Cluster is false

      int formedClusters = -1;

      const int cluSizeIn = clus.size();

      if (digitMCTruth == nullptr) {
        LOGP(info, " Cluster formed Check not using MCTruth");
        formedClusters = clu.solve(&clus, pUserCut, isUnfold); // solve this cluster and add all unfolded clusters to provided list
      }

      // for setting cluster MC-labels
      if (digitMCTruth != nullptr) {
        LOGP(info, " Cluster formed Check mcTruth");

        std::map<int, std::vector<int>> resolvedIndicesMap;

        // solving convoluted clusters, and getting local indices
        formedClusters = clu.solveMC(&clus, pUserCut, isUnfold, resolvedIndicesMap); // solve this cluster and add all unfolded clusters to provided list

        const int startIndex = cluSizeIn;

        const int cluSizeOut = clus.size();
        const int endIndex = cluSizeOut;
        LOGP(info, "clu solve ok : cluSizeIn {} -> cluSizeOut {} formedClusters {}", cluSizeIn, cluSizeOut, formedClusters);
        LOGP(info, "Looping over resolved clusters");
        LOGP(info, "startIndex {} -> endIndex {} formedClusters {}", startIndex, endIndex, formedClusters);
        if (formedClusters > 1 && formedClusters < 6) {
          // use resolvedIndicesMap
          // for(int i = 0; i < resolvedClusters; i++) {
          //  cosnt auto& cluster = clus[i+cluSizeIn];
          //  vector<int> resolvedIndices = resolvedIndicesMap[i];
          //

          // map from "local" to global

   
   				LOGP(info, "\n\n\n=======================\n=======================");
          for (const auto& resolvedIndices : resolvedIndicesMap) {
            // 0.1.2....
            int numSolvedClu = resolvedIndices.first;

            // get cluster
            LOGP(info, "\n=======================\n=======================");
            
            LOGP(info, "Resolved cluster {} of {}", numSolvedClu + 1, formedClusters);

            int currentEntry = startIndex + numSolvedClu;
            LOGP(info, "currentEntry {} startIndex {} ", currentEntry, startIndex);

            const auto& cluster = clus[currentEntry];

            std::vector<int> resolvedInds = resolvedIndices.second;
            // get all the digits, and all the MC-labels from them

            LOGP(info, "clu {}/{}",resolvedIndices.first, clus.size());
	    // iterate over the local to map to global and fidn selected
            
            std::vector<int> globalInd;            
            for (const auto& index : resolvedIndices.second) {
              globalInd.push_back(indicesUnresolved[index]);
            }
            
            // ef: add catch if globalInd.size == 0 ? and add empty label?
            if(globalInd.size() == 0) {
              LOGP(warn, "WARNING globalInd had no entries!");
            }
            
						const int nEntrisIn = mClsLabels->getIndexedSize();
		        iterateMcEntries(cluster, digs, globalInd, digitMCTruth, mClsLabels, clus.size());
		        
		        const int nEntriesOut = mClsLabels->getIndexedSize();
		        
		        LOGP(info, "LOOP number of clusters {} | ndigsResolved {}", clus.size(), cluster.size());
		        
		        LOGP(info, "LOOP IN {} OUT {} ; number of objs in cluster-headArray ", nEntrisIn, nEntriesOut);                      
                                            
          }
					
          if(mClsLabels) {
            const int nEntriesOut = static_cast<int>(mClsLabels->getIndexedSize());
		  			int diff = nEntriesOut - clus.size();

		  			
		  			if(diff!=0) {
		  				Printf("ulik diff %d", diff);

				      if (nEntriesOut != clu.size()) {
				  			LOGP(info, "formedClusters = {} unequal sizes of labels ({}) and clusters ({}) diff : {}", formedClusters, nEntriesOut, clus.size(), diff);
				  		 Printf("diff %d", diff);
							}

				      if (nEntriesOut == clu.size()) {
				  			LOGP(info, "formedClusters = {} equal sizes of labels ({}) and clusters ({}) diff : {}", formedClusters, nEntriesOut, clus.size(), diff);
							}
		  			}	
		  			
		  			else {
		  				Printf("diff ok %d", diff);
		  			}	  					  				  		
    			}
        } else if (formedClusters == 1 || formedClusters >= 6 ) { // 6 is max allowed value kMaxLocMax
          LOGP(info, "\n\n\n=======================\n =======================");
          //  const auto& cluster = clus[i+cluSizeIn];
          // the "unresolved" are here also the resolved,
          // since we only had 1 cluster to begin with in raw cluster
          // we can directly use indicesUnresolved as the global
          const auto& cluster = clus[startIndex];
          /*for (const auto& index : indicesUnresolved) {
              const auto& dig = digs[index];
              // index is used for mcDigitsTruth->getLabels(index)
              gsl::span<const o2::MCCompLabel> mcArray = digitMCTruth->getLabels(index);
          }*/
					const int nEntrisIn = mClsLabels->getIndexedSize();
          iterateMcEntries(cluster, digs, indicesUnresolved, digitMCTruth, mClsLabels, clus.size());
					
					LOGP(info, "number of clusters {} | ndigsResolved {}", clus.size(), cluster.size());
          
					
          if(mClsLabels) {
             const int nEntriesOut = static_cast<int>(mClsLabels->getIndexedSize());
					LOGP(info, "IN {} OUT {} ; number of objs in cluster-headArray ", nEntrisIn, nEntriesOut);
					
		  			int diff = nEntriesOut - clus.size();
		  			
		  			if(diff!=0) {
		  				Printf("ulik diff %d", diff);
				      if (nEntriesOut != clu.size()) {
				  			LOGP(info, "formedClusters = {} unequal sizes of labels ({}) and clusters ({}) diff : {}", formedClusters, nEntriesOut, clus.size(), diff);
				  		 Printf("diff %d", diff);
							}

				      if (nEntriesOut == clu.size()) {
				  			LOGP(info, "formedClusters = {} equal sizes of labels ({}) and clusters ({}) diff : {}", formedClusters, nEntriesOut, clus.size(), diff);
							}	
		  			}	else {
		  				Printf("diff ok %d", diff);
		  			}	  				  			
    			}
          
        } else {
          LOGP(info, "\n\n\n formedClusters? {}", formedClusters);
        }
        
        
      } // if digitsMcTruth




    } // digits loop for current chamber

    vPad.clear();

  } // chambers loop

  // const int numCluStart = clus.size();

  LOGP(info, "\n\n ============================== called Dig2Clu\n looping over digits");

  LOGP(info, " new Clusters  {}", clus.size() - numCluStart);

  LOGP(info, " digs Size {}", digs.size());

  return;

} // Dig2Clu()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Clusterer::iterateMcEntries(const Cluster& cluster, gsl::span<const o2::hmpid::Digit> digits, const std::vector<int>& indices, MCLabelContainer const* digitMCTruth, MCLabelContainer* mClsLabels, int cluSize)
{


  LOGP(info, "Based on pos : pdg {} mother {} tid {}", cluster.getPDG(), cluster.getMotherId(), cluster.getTrackId());

  int lbl = mClsLabels->getIndexedSize(); // this should correspond to the current number of clusters? ;

  LOGP(info, "lbl = {} indices_size {} : cluster.size {}", lbl, indices.size(), cluster.size());

  // LOGP(info, "number of clusters {}", clus.size());

  //LOGP(info, "number of objs in cluster-headArray : {}", mClsLabels->getIndexedSize());

  //LOGP(info, "number of objs in cluster-TruthArray : {}", mClsLabels->getNElements());
  // loop over resolvedIndexes array;
  int numDigitsCount = 0;
  for (int i = 0; i < indices.size(); i++) {
    const int& indexOfDigGlobal = indices[i];
    const auto& digOfClu = digits[indexOfDigGlobal];
    LOGP(info, "digit  {}/{} : indexOfDigGlobal {}", i + 1, indices.size(), indexOfDigGlobal);
    

    LOGP(info, "digit numGlobal {} : x {} y {}", indexOfDigGlobal, digOfClu.getX(), digOfClu.getY());
    /*if (iDig < 0 || iDig >= digs.size()) {
      LOGP(info, "iDig out of bounds");
      break;
    }*/
    numDigitsCount++;
    // const auto& digOfClu = &digs[iDig];

    if (digitMCTruth == nullptr) {
      LOGP(info, "digitMCTruth was nullptr");
      continue;
    }

    // ef: TODO: is this correct?
    int digitLabel = digOfClu.getLabel();
    // ef: TODO it should be the index of the digit in the array of digits? Like this :

    // for MC-truth
    digitLabel = indexOfDigGlobal;

    const int digEventNum = digOfClu.getEventNumber();

    const int pdgOfDig = digOfClu.getPDG();

    // printf("digitLabel = %d\n", digitLabel);

    // all MC-hits from the specific digit

    gsl::span<const o2::MCCompLabel> mcArray = digitMCTruth->getLabels(digitLabel);

    // ef : vi maa sette en logikk som setter mcTruth fra digit til riktig cluster-indeks i mClsLabels/clus

    //LOGP(info, "contributing digit = ({}/{}), digitLabel  = {} || pdg of digit {}", numDigitsCount, cluster.size(), digitLabel, pdgOfDig);

    LOGP(info, "mcArray size {}", mcArray.size());
    LOGP(info, "======= Looping Labels of dig ===== ");
    for (int j = 0; j < static_cast<int>(mcArray.size()); j++)

    {

      const auto& currentIndex = digitMCTruth->getMCTruthHeader(digitLabel).index + j;
      auto label = digitMCTruth->getElement(currentIndex);

      //LOGP(info, "digitLabel {}, digitMCTruth->getMCTruthHeader(digitLabel).index {}, j {}", digitLabel, digitMCTruth->getMCTruthHeader(digitLabel).index, j);
      //LOGP(info, "digitMCTruth->getElement({}/{})", currentIndex, digitMCTruth->getIndexedSize());

      // same as digits having multiple hits
      // clsuters have multiple digits
      // we fill MC-Complabel label at index lbl for headArray
      // this is the hit for a digit in the cluster

      LOGP(info, "adding mc label at index {}", lbl);

      // skal ikke dette være lblFromClu??

      // lbl = mClsLabels->getIndexedSize()

      auto lblFromClu = cluSize - 1;

      // Det ga runtime error

      // må vi gjøre dette i clu. solve?

      // lbl, label, mcArray

      mClsLabels->addElement(lbl, label);
      LOGP(info, "number of labels for clu {}  : lbl", mClsLabels->getLabels(lbl).size());
      /*
      LOGP(info, "number of clusters {}", cluSize);

      LOGP(info, "number of objs in cluster-headArray : {}", mClsLabels->getIndexedSize());

      LOGP(info, "number of objs in cluster-TruthArray : {}", mClsLabels->getNElements());

      LOGP(info, "number of labels for clu {} : clus.size() - 1", mClsLabels->getLabels(cluSize - 1).size());

      */ 

      /*
      /// query an MC track given a basic label object
      /// returns nullptr if no track was found
      MCTrack const* getTrack(o2::MCCompLabel const&) const;
      */
      // ef : this means we get the track of the hit
      // ef :TODO remove print statements or add if

      const o2::MCTrack* mcTrack = nullptr;

      const o2::MCTrack* mcTrackFromDig = nullptr;

      const o2::MCTrack* mcTrackFromMother = nullptr;

      const auto& mcReader = std::make_unique<o2::steer::MCKinematicsReader>("collisioncontext.root");

      bool printVals = true;

      if (printVals) {
        if (!mcReader)

        {

          LOGP(info, "mcReader nullptr");

          continue;
        }

        if (mcReader->getTrack(label))

        {

          try

          {

            mcTrack = mcReader->getTrack(label);
          }

          catch (const std::exception& e)

          {

            LOGP(error, "       Exception caught while trying to read MC track: %s", e.what());

            continue;
          }

          catch (...)

          {

            LOGP(error, "       Unknown exception caught while trying to read MC track");

            continue;
          }
        }

        else

        {

          LOGP(info, "mcReader->getTrack(label) gave nullptr");
        }

        int pdgDigMcTruth = -2, pdgDigit = -2, pdgMother = -2;

        try

        {

          pdgDigMcTruth = mcTrack->GetPdgCode();
        }

        catch (const std::exception& e)

        {

          LOGP(error, "       Exception caught while trying to read MC track: %s", e.what());
        }

        catch (...)

        {

          LOGP(error, "       Unknown exception caught while trying to read MC track");
        }

        // TParticlePDG* pPDG = TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdgCode());

        int trackID, evID, srcID;

        bool fake;

        const auto eid = digOfClu.getEventNumber();

        const auto tid = digOfClu.getTrackId();

        const auto mid = digOfClu.getMotherId();

        const auto sid = digOfClu.getSourceId();

        if (!mcReader->getTrack(eid, tid))

        {

          LOGP(info, "nullptr mcReader->getTrack(eid, tid)");

          try

          {

            mcTrackFromDig = mcReader->getTrack(eid, tid);

            pdgDigit = mcTrackFromDig->GetPdgCode();
          }

          catch (const std::exception& e)

          {

            LOGP(error, "       Exception caught while trying to read MC track: %s", e.what());
          }

          catch (...)

          {

            LOGP(error, "       Unknown exception caught while trying to read MC track");
          }
        }

        else

        {

          LOGP(debug, "mcReader->getTrack(eid, tid) gave nullptr");
        }

        if (!mcReader->getTrack(eid, mid))

        {

          // LOGP(info, "nullptr mcReader->getTrack(eid, mid)");

          try

          {

            mcTrackFromDig = mcReader->getTrack(eid, mid);

            pdgMother = mcTrackFromMother->GetPdgCode();
          }

          catch (const std::exception& e)

          {

            LOGP(error, "       Exception caught while trying to read MC track: %s", e.what());
          }

          catch (...)

          {

            LOGP(error, "       Unknown exception caught while trying to read MC track");
          }
        }

        trackID = mcArray[j].getTrackID(); // const { return static_cast<int>(mLabel & maskTrackID); }

        evID = mcArray[j].getEventID(); // const { return isFake() ? -getTrackID() : getTrackID(); }

        srcID = mcArray[j].getSourceID(); // const { return (mLabel >> nbitsTrackID) & maskEvID; }

        fake = mcArray[j].isFake(); // const { return (mLabel >> (nbitsTrackID + nbitsEvID)) & maskSrcID; }

        // ef : nt marked as const in MCOMPLabel

        /// mcArray[j].get(trackID, evID, srcID, fake);

        // LOGP(info, "checking element {}/{} in array of labels", j + 1, mcArray.size());

        

        LOGP(info, "mcArray : evID {}, trackID {}, mid {}, srcID {}, fake {}", evID, trackID, mcTrack->getMotherTrackId(), srcID, fake);

        LOGP(info, "from digit : eid {}, tid {}, mid {}, sid {}", eid, tid, mid, sid);
				/*
        // printf("checking element %d in the array of labels\n", j);

        LOGP(info, "EventID from MC-label = {}; from dig : {}", evID, digEventNum);

        */

        LOGP(info, " pdg of digit {} || MC digit label {} | pdg from (digit-eid, digit-tid) {} | of mother (digit-eid, digit-mid) {}", pdgOfDig, pdgDigMcTruth, pdgDigit, pdgMother);
	if(pdgOfDig!=pdgDigMcTruth) {
	  LOGP(info, "pdgOfDig ulik pdgDigMcTruth");
	}
        // LOGP(info, "num Digits : = {}", digs.size());

      } // end if printVals
    }   // end for mcArray
  }     // end for iDig : resolvedIndices
  /* 
  if (mClsLabels) {
    LOGP(info, "number of objs in cluster-headArray : {}", mClsLabels->getIndexedSize());
    LOGP(info, "number of objs in cluster-TruthArray : {}", mClsLabels->getNElements());
    
  } */ // end if mClsLabels

} // end iterateMcEntries
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Clusterer::FormClu(Cluster& pClu, int pDig, gsl::span<const o2::hmpid::Digit> digs, TMatrixF& pDigMap)

{

  // Forms the initial cluster as a combination of all adjascent digits. Starts from the given digit then calls itself recursevly  for all neighbours.

  // Arguments: pClu - pointer to cluster being formed

  //   Returns: none

  //   void digAdd(const o2::hmpid::Digit* pDig);                                         // add new digit to the cluster

  pClu.digAdd(&digs[pDig]); // take this digit in cluster

  // ef : add index of digit to cluster

  // pClu.setUnresolvedIndex(pDig);

  LOGP(info, "FormClu : pDig = {} | x {} y {}", pDig, digs[pDig].getX(), digs[pDig].getY());

  int cnt = 0;

  int cx[4];

  int cy[4];

  int padChX = 0;

  int padChY = 0;

  int module = 0;

  o2::hmpid::Digit::pad2Absolute(digs[pDig].getPadID(), &module, &padChX, &padChY);

  if (padChX > Param::kMinPx)

  {

    cx[cnt] = padChX - 1;

    cy[cnt] = padChY;

    cnt++;

  } // left

  if (padChX < Param::kMaxPcx)

  {

    cx[cnt] = padChX + 1;

    cy[cnt] = padChY;

    cnt++;

  } // right

  if (padChY > Param::kMinPy)

  {

    cx[cnt] = padChX;

    cy[cnt] = padChY - 1;

    cnt++;

  } // down

  if (padChY < Param::kMaxPcy)

  {

    cx[cnt] = padChX;

    cy[cnt] = padChY + 1;

    cnt++;

  } // up

  for (int i = 0; i < cnt; i++)

  { // neighbours loop

    pDig = UseDig(cx[i], cy[i], pDigMap);

    if (pDig != -1)

    {

      FormClu(pClu, pDig, digs, pDigMap);

    } // check if this neighbour pad fired and mark it as taken

  } // neighbours loop

} // FormClu()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Clusterer::FormCluMC(Cluster& pClu, int pDig, gsl::span<const o2::hmpid::Digit> digs, TMatrixF& pDigMap, std::vector<int>& indicesUnresolved)

{

  // Forms the initial cluster as a combination of all adjascent digits. Starts from the given digit then calls itself recursevly  for all neighbours.

  // Arguments: pClu - pointer to cluster being formed

  //   Returns: none

  //   void digAdd(const o2::hmpid::Digit* pDig);                                         // add new digit to the cluster

  pClu.digAdd(&digs[pDig]); // take this digit in cluster
  indicesUnresolved.push_back(pDig);
  // ef : add index of digit to cluster

  // pClu.setUnresolvedIndex(pDig);

  LOGP(info, "FormClu : pDig = {} | x {} y {}", pDig, digs[pDig].getX(), digs[pDig].getY());

  int cnt = 0;

  int cx[4];

  int cy[4];

  int padChX = 0;

  int padChY = 0;

  int module = 0;

  o2::hmpid::Digit::pad2Absolute(digs[pDig].getPadID(), &module, &padChX, &padChY);

  if (padChX > Param::kMinPx)

  {

    cx[cnt] = padChX - 1;

    cy[cnt] = padChY;

    cnt++;

  } // left

  if (padChX < Param::kMaxPcx)

  {

    cx[cnt] = padChX + 1;

    cy[cnt] = padChY;

    cnt++;

  } // right

  if (padChY > Param::kMinPy)

  {

    cx[cnt] = padChX;

    cy[cnt] = padChY - 1;

    cnt++;

  } // down

  if (padChY < Param::kMaxPcy)

  {

    cx[cnt] = padChX;

    cy[cnt] = padChY + 1;

    cnt++;

  } // up

  for (int i = 0; i < cnt; i++)

  { // neighbours loop

    pDig = UseDig(cx[i], cy[i], pDigMap);

    if (pDig != -1)

    {

      FormCluMC(pClu, pDig, digs, pDigMap, indicesUnresolved);

    } // check if this neighbour pad fired and mark it as taken

  } // neighbours loop

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

  if (pDig->getQ() > 4.)

  {

    return true;
  }

  else

  {

    return false;
  }
}


