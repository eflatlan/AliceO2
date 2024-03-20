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

      LOGP(info, "\n\n ================================\n New iDig : {} \n================================\n\n", iDig);

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

      // ef : ad digitMCTruth

      FormClu(clu, pUsedDig, digs, padMap); // form cluster starting from this digit by recursion
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

      // clu.solveMC(&clus, pUserCut, isUnfold, mClsLabels, digitMCTruth);
      // int formedClusters = -1;
      int formedClusters = clu.solve(&clus, pUserCut, isUnfold); // solve this cluster and add all unfolded clusters to provided list

      LOGP(info, " Cluster formed Check mcTruth");

 

      // ,  MCLabelContainer const* digitMCTruth

      // denne skal ikke vaere nptr

      /*

      if(clu.dig(0)==nullptr)

      Printf("dig2Clu Digit i0 nullptr ");


      else


      Printf("dig2Clu Digit i0 ok");

      */

      LOGP(info, "set CluCh {}", iCh);

      clu.setCh(iCh);

      // ef: TODO remvoe

      const int cluSizeIn = clus.size();


      const int startIndex = cluSizeIn;
      const int cluSizeOut = clus.size();




      // ef : TODO: move MCLabels setting here ?

      if (digitMCTruth != nullptr) {
        // ef : TODO: move MCLabels setting here ?
        // clu.setMCTruth(digitMCTruth); // we do this in Cluster.h

        // just set it normmaly
        /*if(formedClusters == 1) {

        }* /
        // we need to set the correct digits for the clusters
        // and the correct MC-labels from them
        else if (formedClusters > 1) {



          for(int i = startIndex; i < cluSizeOut; i++)
          {
            auto& cluster = clus[i];

            // if useMC : clean pointers after cluster is added to loop
            cluster.cleanPointer()
            // clusters are cleaned here, so we do not have acces to digit-pointers?
          }*/

        const int cluSizeOut = clus.size();

        LOGP(info, "clu solve ok : cluSizeIn {} -> cluSizeOut {} formedClusters{}", cluSizeIn, cluSizeOut, formedClusters);

        const int endIndex = cluSizeOut;

        // Loop over resolved clusters
        // if only 1, it will still loop once
        LOGP(info, "Looping over resolved clusters");
        LOGP(info, "startIndex {} -> endIndex {} formedClusters{}", startIndex, endIndex, formedCl);
        int resolvdCluCnt = 0;
        for (int i = startIndex; i < endIndex; i++) {
          LOGP(info, "Resolved cluster {} of {}", resolvdCluCnt++, formedClustes);

          const auto& cluster = clus[i];
          //const auto& resovledIndexes = getResolvedIndexes();
          LOGP(info, "index {}, cluster size {})", i, cluster.size());

          LOGP(info, "Based on pos : pdg {} mother {} tid {}", cluster.getPDG(), cluster.getMotherId(), cluster.getTrackId());


          int lbl = mClsLabels->getIndexedSize(); // this should correspond to the current number of clusters? ;

          LOGP(info, "lbl = {} (mClsLabels->getIndexedSize())", lbl);

          LOGP(info, "number of clusters {}", clus.size());

          LOGP(info, "number of objs in cluster-headArray : {}", mClsLabels->getIndexedSize());

          LOGP(info, "number of objs in cluster-TruthArray : {}", mClsLabels->getNElements());
          
          
          // loop over resovledIndexes array; 
          int numDigitsCount = 0;
          int indexOfFirst = resovledIndexes[0];
          
          // special case: we have only one digit in the cluster
          // and it is in the first position of the digit-array
          /*if(indexOfFirst == 0 && cluster.size() == 1) {          
          } */        


          // resovledIndexes is an array containing the indices of the
          // digits from digs that were used to form the cluster
          for(const auto& iDig : resovledIndexes) {
            
            // iDig < -1 : means we have empty digit
            if (iDig < 0 || iDig >= digs.size()) {
              LOGP(info, "iDig out of bounds");
              break;
            }

            const auto& digOfClu = digs[iDig];

            if(!digOfClu) {
              LOGP(info, "digOfClu was nullptr");
              continue;
            }

            const auto& digMCTruth = digOfClu->getMCTruth();
            if(digMCTruth == nullptr) {
              LOGP(info, "digMCTruth was nullptr");
              continue;
            }  
            
            // ef: TODO: is this correct?

            int digitLabel = digOfClu->getLabel();

            // ef: TODO it should be the index of the digit in the array of digits? Like this :

            digitLabel = iDig;

            const int digEventNum = dig->getEventNumber();

            const int pdgOfDig = digOfClu->getPDG();

            // printf("digitLabel = %d\n", digitLabel);

            // all MC-hits from the specific digit
            gsl::span<const o2::MCCompLabel> mcArray = digitMCTruth->getLabels(digitLabel);
            // ef : vi maa sette en logikk som setter mcTruth fra digit til riktig cluster-indeks i mClsLabels/clus

            LOGP(info, "contributing digit = ({}/{}), digitLabel  = {} || pdg of digit {}", digIndex + 1, digsOfClu.size(), digitLabel, pdgOfDig);

            LOGP(info, "mcArray size {}", mcArray.size());

            LOGP(info, "======= Looping Labels of dig ===== ");

            for (int j = 0; j < static_cast<int>(mcArray.size()); j++)

            {

              const auto& currentIndex = digitMCTruth->getMCTruthHeader(digitLabel).index + j;

              auto label = digitMCTruth->getElement(currentIndex);

              LOGP(info, "digitLabel {}, digitMCTruth->getMCTruthHeader(digitLabel).index {}, j {}", digitLabel, digitMCTruth->getMCTruthHeader(digitLabel).index, j);

              LOGP(info, "digitMCTruth->getElement({}/{})", currentIndex, digitMCTruth->getIndexedSize());

              // same as digits having multiple hits

              // clsuters have multiple digits

              // we fill MC-Complabel label at index lbl for headArray

              // this is the hit for a digit in the cluster

              LOGP(info, "adding mc label at index {}", lbl);

              // skal ikke dette være lblFromClu??

              // lbl = mClsLabels->getIndexedSize()

              auto lblFromClu = clus.size() - 1;

              // Det ga runtime error

              // må vi gjøre dette i clu. solve?

              // lbl, label, mcArray

              mClsLabels->addElement(lbl, label);

              LOGP(info, "number of clusters {}", clus.size());

              LOGP(info, "number of objs in cluster-headArray : {}", mClsLabels->getIndexedSize());

              LOGP(info, "number of objs in cluster-TruthArray : {}", mClsLabels->getNElements());

              LOGP(info, "number of labels for clu {} : clus.size() - 1", mClsLabels->getLabels(clus.size() - 1).size());

              LOGP(info, "number of labels for clu {}  : lbl", mClsLabels->getLabels(lbl).size());

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

              if (printVals)

              {

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

                const auto eid = digOfClu->getEventNumber();

                const auto tid = digOfClu->getTrackId();

                const auto mid = digOfClu->getMotherId();

                const auto sid = digOfClu->getSourceId();

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

                  LOGP(info, "mcReader->getTrack(eid, tid) gave nullptr");
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

                LOGP(info, "checking element {}/{} in array of labels", j + 1, mcArray.size());

                /*

                LOGP(info, "mcArray : evID {}, trackID {}, srcID {}, fake {}", evID, trackID, srcID, fake);



                LOGP(info, "from digit : eid {}, tid {}, mid {}, sid {}", eid, tid, mid, sid);



                // printf("checking element %d in the array of labels\n", j);



                LOGP(info, "EventID from MC-label = {}; from dig : {}", evID, digEventNum);

                */

                LOGP(info, " pdg of digit {} || MC digit label {} | pdg from (digit-eid, digit-tid) {} | of mother (digit-eid, digit-mid) {}", pdgOfDig, pdgDigMcTruth, pdgDigit, pdgMother);

                // LOGP(info, "num Digits : = {}", digs.size());

              } // end if printVals

            } // end for mcArray

          } // end for digIndex
          // if useMC : clean pointers after cluster is added to loop
          // cluster.cleanPointer()
          // clusters are cleaned here, so we do not have acces to digit-pointers?
        } // end resolved clusters-loop

        if (mClsLabels) {

          LOGP(info, "number of objs in cluster-headArray : {}", mClsLabels->getIndexedSize());

          LOGP(info, "number of objs in cluster-TruthArray : {}", mClsLabels->getNElements());

          if (mClsLabels->getIndexedSize() != cluSizeOut) {

            LOGP(info, "unequal sizes of labels and clusters");
          }
        } // end if mClsLabels


      } // end if digitMCTruth
      /*

      if(clus.back().dig(0) == nullptr) {Printf("dig2Clu dig was nullptr!!");}



      */

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

void Clusterer::FormClu(Cluster& pClu, int pDig, gsl::span<const o2::hmpid::Digit> digs, TMatrixF& pDigMap)

{

  // Forms the initial cluster as a combination of all adjascent digits. Starts from the given digit then calls itself recursevly  for all neighbours.

  // Arguments: pClu - pointer to cluster being formed

  //   Returns: none

  //   void digAdd(const o2::hmpid::Digit* pDig);                                         // add new digit to the cluster

  pClu.digAdd(&digs[pDig]); // take this digit in cluster

  // ef : add index of digit to cluster
  pClu.setUnresolvedIndex(pDig);
  LOGP(info, "FormClu : pDig = {} | x {} y {}", pDig, digs[pDig].x(), digs[pDig].y());

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
