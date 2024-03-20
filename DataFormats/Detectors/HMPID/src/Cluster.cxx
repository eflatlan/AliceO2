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

#include <TRandom.h>

#include <TMarker.h>

#include "DataFormatsHMP/Cluster.h"

#include "Framework/Logger.h"

#include <TGeoManager.h>

#include <TVirtualFitter.h>

#include <cmath>

ClassImp(o2::hmpid::Cluster);

namespace o2

{

namespace hmpid

{

bool o2::hmpid::Cluster::fgDoCorrSin = true;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Calculates naive cluster position as a center of gravity of its digits.

void Cluster::coG()

{

  if (!mDigs) {

    LOGP(fatal, "digits are missing in the cluster");
  }

  int minPadX = 999;

  int minPadY = 999;

  int maxPadX = -1;

  int maxPadY = -1; // for box finding

  if (mDigs->size() == 0) {

    return;

  } // no digits in this cluster

  mXX = mYY = mQRaw = 0; // init summable parameters

  mCh = -1; // init chamber

  int maxQpad = -1;

  int maxQ = -1; // to calculate the pad with the highest charge

  o2::hmpid::Digit* pDig = nullptr;

  for (int iDig = 0; iDig < mDigs->size(); iDig++) { // digits loop

    int x, y, mod;

    int padId = (*mDigs)[iDig]->getPadID();

    o2::hmpid::Digit::pad2Absolute(padId, &mod, &x, &y);

    if (x > maxPadX) {

      maxPadX = x;

    } // find the minimum box that contain the cluster  MaxX

    if (y > maxPadY) {

      maxPadY = y;

    } // MaxY

    if (x < minPadX) {

      minPadX = x;

    } // MinX

    if (y < minPadY) {

      minPadY = y;

    } // MinY

    float q = (*mDigs)[iDig]->mQ; // get QDC

    mXX += o2::hmpid::Digit::lorsX(padId) * q;

    mYY += o2::hmpid::Digit::lorsY(padId) * q; // add digit center weighted by QDC

    mQRaw += q; // increment total charge

    if (q > maxQ) {

      maxQpad = padId;

      maxQ = (int)q;

    } // to find pad with highest charge

    mCh = mod; // initialize chamber number

  } // digits loop

  mBox = (maxPadX - minPadX + 1) * 100 + maxPadY - minPadY + 1; // dimension of the box: format Xdim*100+Ydim

  if (mQRaw != 0) {

    mXX /= mQRaw;

    mYY /= mQRaw;

  } // final center of gravity

  if (mDigs->size() > 1 && fgDoCorrSin) {

    corrSin();

  } // correct it by sinoid

  mQ = mQRaw; // Before starting fit procedure, Q and QRaw must be equal

  mMaxQpad = maxQpad;

  mMaxQ = maxQ; // store max charge pad to the field

  mChi2 = 0; // no Chi2 to find

  mNlocMax = 0; // proper status from this method

  mSt = o2::hmpid::Cluster::kCoG;

  return;

  //  setClusterParams(mXX, mYY, mCh); //need to fill the AliCluster3D part

} // CoG()

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Correction of cluster x position due to sinoid, see HMPID TDR  page 30

void Cluster::corrSin()

{

  const auto param = o2::hmpid::Param::instanceNoGeo();

  int pc;

  int px;

  int py;

  o2::hmpid::Param::lors2Pad(mXX, mYY, pc, px, py); // tmp digit to get it center

  double x = mXX - param->lorsX(pc, px); // diff between cluster x and center of the pad contaning this cluster

  double xpi8on10 = M_PI / 0.8 * x;

  mXX += 3.31267e-2 * sin(2.0 * xpi8on10) - 2.66575e-3 * sin(4.0 * xpi8on10) + 2.80553e-3 * sin(6.0 * xpi8on10) + 0.0070;

  return;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*void Cluster::draw(Option_t*)

{

  TMarker *pMark = new TMarker(X(), Y(), 5);

  pMark->SetUniqueID(mSt);pMark->SetMarkerColor(kBlue); pMark->Draw();

}

*/

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Cluster::fitFunc(int& iNpars, double* deriv, double& chi2, double* par, int iflag)

{

  // Cluster fit function

  // par[0]=x par[1]=y par[2]=q for the first Mathieson shape

  // par[3]=x par[4]=y par[5]=q for the second Mathieson shape and so on up to iNpars/3 Mathieson shapes

  // For each pad of the cluster calculates the difference between actual pad charge and the charge induced to this pad by all Mathieson distributions

  // Then the chi2 is calculated as the sum of this value squared for all pad in the cluster.

  // Arguments: iNpars - number of parameters which is number of local maxima of cluster * 3

  //            chi2   - function result to be minimised

  //            par   - parameters array of size iNpars

  //   Returns: none

  Cluster* pClu = (Cluster*)TVirtualFitter::GetFitter()->GetObjectFit();

  int nPads = pClu->mSi;

  chi2 = 0;

  int iNshape = iNpars / 3;

  for (int i = 0; i < nPads; i++) { // loop on all pads of the cluster

    double dQpadMath = 0;

    for (int j = 0; j < iNshape; j++) { // Mathiesons loop as all of them may contribute to this pad

      int baseOff = 3 * j;

      int baseOff1 = baseOff + 1;

      int baseOff2 = baseOff + 2;

      double fracMathi = o2::hmpid::Digit::intMathieson(par[baseOff], par[baseOff1], pClu->dig(i)->getPadID());

      dQpadMath += par[baseOff2] * fracMathi; // par[3*j+2] is charge par[3*j] is x par[3*j+1] is y of current Mathieson
    }

    if (dQpadMath > 0 && pClu->dig(i)->mQ > 0) {

      chi2 += std::pow((pClu->dig(i)->mQ - dQpadMath), 2.0) / pClu->dig(i)->mQ; // chi2 function to be minimized
    }
  }

  //---calculate gradients...

  if (iflag == 2) {

    std::vector<std::vector<float>> derivPart(iNpars, std::vector<float>(nPads, 0.0f));

    double dHalfPadX = o2::hmpid::Param::sizeHalfPadX();

    double dHalfPadY = o2::hmpid::Param::sizeHalfPadY();

    for (int i = 0; i < nPads; i++) { // loop on all pads of the cluster

      int iPadId = pClu->dig(i)->getPadID();

      double lx = o2::hmpid::Digit::lorsX(iPadId);

      double ly = o2::hmpid::Digit::lorsY(iPadId);

      for (int j = 0; j < iNshape; j++) { // Mathiesons loop as all of them may contribute to this pad

        int baseOff = 3 * j;

        int baseOff1 = baseOff + 1;

        int baseOff2 = baseOff + 2;

        double fracMathi = o2::hmpid::Digit::intMathieson(par[baseOff], par[baseOff1], iPadId);

        derivPart[baseOff][i] += par[baseOff2] * (o2::hmpid::Digit::mathiesonX(par[baseOff] - lx - dHalfPadX) - o2::hmpid::Digit::mathiesonX(par[baseOff] - lx + dHalfPadX)) *

                                 o2::hmpid::Digit::intPartMathiY(par[baseOff1], iPadId);

        derivPart[baseOff1][i] += par[baseOff2] * (o2::hmpid::Digit::mathiesonY(par[baseOff1] - ly - dHalfPadY) - o2::hmpid::Digit::mathiesonY(par[baseOff1] - ly + dHalfPadY)) *

                                  o2::hmpid::Digit::intPartMathiX(par[baseOff], iPadId);

        derivPart[baseOff2][i] += fracMathi;
      }
    }

    // loop on all pads of the cluster

    for (int i = 0; i < nPads; i++) { // loop on all pads of the cluster

      int iPadId = pClu->dig(i)->getPadID();

      double dPadmQ = pClu->dig(i)->mQ;

      double dQpadMath = 0.0; // pad charge collector

      double twoOverMq = 2.0 / dPadmQ;

      for (int j = 0; j < iNshape; j++) { // Mathiesons loop as all of them may contribute to this pad

        int baseOff = 3 * j;

        int baseOff1 = baseOff + 1;

        int baseOff2 = baseOff + 2;

        double fracMathi = o2::hmpid::Digit::intMathieson(par[baseOff], par[baseOff1], iPadId);

        dQpadMath += par[baseOff2] * fracMathi;

        if (dQpadMath > 0 && dPadmQ > 0) {

          double appoggio = twoOverMq * (dPadmQ - dQpadMath);

          deriv[baseOff] += appoggio * derivPart[baseOff][i];

          deriv[baseOff1] += appoggio * derivPart[baseOff1][i];

          deriv[baseOff2] += appoggio * derivPart[baseOff2][i];
        }
      }
    }

  } //---gradient calculations ended

  // fit ended. Final calculations

} // FitFunction()

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Print current cluster

void Cluster::print(Option_t* opt) const

{

  const char* status = nullptr;

  switch (mSt) {

    case kFrm:

      status = "formed        ";

      break;

    case kUnf:

      status = "unfolded (fit)";

      break;

    case kCoG:

      status = "coged         ";

      break;

    case kLo1:

      status = "locmax 1 (fit)";

      break;

    case kMax:

      status = "exceeded (cog)";

      break;

    case kNot:

      status = "not done (cog)";

      break;

    case kEmp:

      status = "empty         ";

      break;

    case kEdg:

      status = "edge     (fit)";

      break;

    case kSi1:

      status = "size 1   (cog)";

      break;

    case kNoLoc:

      status = "no LocMax(fit)";

      break;

    case kAbn:

      status = "Abnormal fit  ";

      break;

    case kBig:

      status = "Big Clu(>100) ";

      break;

    default:

      status = "??????";

      break;
  }

  float ratio = 0;

  if (mQ > 0.0 && mQRaw > 0.0) {

    ratio = mQ / mQRaw * 100;
  }

  printf("%sCLU: ch=%i  (X = %7.3f, Y = %7.3f) Q=%8.3f Qraw=%8.3f(%3.0f%%) Size=%2i DimBox=%i LocMax=%i Chi2=%7.3f   %s\n",

         opt, mCh, mXX, mYY, mQ, mQRaw, ratio, mSi, mBox, mNlocMax, mChi2, status);

  if (mDigs && mDigs->size() > 0) {

    std::cout << "Digits of Cluster" << std::endl;

    for (int i; i < mDigs->size(); i++) {

      std::cout << (*mDigs)[i] << std::endl;
    }
  }

  return;

} // Print()

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int Cluster::solve(std::vector<o2::hmpid::Cluster>* pCluLst, float* pSigmaCut, bool isTryUnfold)

{

  // This methode is invoked when the cluster is formed to solve it. Solve the cluster means to try to unfold the cluster

  // into the local maxima number of clusters. This methode is invoked by AliHMPIDRconstructor::Dig2Clu() on cluster by cluster basis.

  // At this point, cluster contains a list of digits, cluster charge and size is precalculated in AddDigit(), position is preset to (-1,-1) in ctor,

  // status is preset to kFormed in AddDigit(), chamber-sector info is preseted to actual values in AddDigit()

  // Method first finds number of local maxima and if it's more then one tries to unfold this cluster into local maxima number of clusters

  // Arguments: pCluLst     - cluster list pointer where to add new cluster(s)

  //            isTryUnfold - flag to switch on/off unfolding

  //   Returns: number of local maxima of original cluster

  const auto param = o2::hmpid::Param::instanceNoGeo();

  if (!mDigs) {

    LOGP(fatal, "digits are missing in the cluster");
  }

  const int kMaxLocMax = 6; // max allowed number of loc max for fitting

  coG(); // First calculate CoG for the given cluster

  int iCluCnt = pCluLst->size(); // get current number of clusters already stored in the list by previous operations

  int rawSize = mSi; // get current raw cluster size

  if (rawSize > 100) {

    mSt = kBig;
  } else if (isTryUnfold == false) {

    mSt = kNot;
  } else if (rawSize == 1) {

    mSt = kSi1;
  }

  if (rawSize > 100 || isTryUnfold == false || rawSize == 1) { // No deconv if: 1 - big cluster (also avoid no zero suppression!)

    // setClusterParams(mXX, mYY, mCh); //                               2 - flag is set to FALSE

    // new ((*pCluLst)[iCluCnt++]) Cluster(*this); //                      3 - size = 1

    pCluLst->push_back(o2::hmpid::Cluster(*this));

    pCluLst->back().cleanPointers();

    return 1; // add this raw cluster
  }

  //  Phase 0. Initialise Fitter

  double arglist[10]{0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  float ierflg = 0.;

  TVirtualFitter* fitter = TVirtualFitter::Fitter((TObject*)this, 3 * 6); // initialize Fitter

  if (fitter == nullptr) {

    LOG(fatal) << "TVirtualFitter could not be created";

    return 1;
  }

  arglist[0] = -1;

  ierflg = fitter->ExecuteCommand("SET PRI", arglist, 1); // no printout

  ierflg = fitter->ExecuteCommand("SET NOW", arglist, 0); // no warning messages

  arglist[0] = 1;

  ierflg = fitter->ExecuteCommand("SET GRA", arglist, 1); // force Fitter to use my gradient

  fitter->SetFCN(Cluster::fitFunc);

  // Phase 1. Find number of local maxima. Strategy is to check if the current pad has QDC more then all neigbours. Also find the box contaning the cluster

  mNlocMax = 0;

  for (int iDig1 = 0; iDig1 < rawSize; iDig1++) { // first digits loop

    auto pDig1 = (*mDigs)[iDig1]; // take next digit

    int iCnt = 0; // counts how many neighbouring pads has QDC more then current one

    for (int iDig2 = 0; iDig2 < rawSize; iDig2++) { // loop on all digits again

      if (iDig1 == iDig2) {

        continue;

      } // the same digit, no need to compare

      auto pDig2 = (*mDigs)[iDig2]; // take second digit to compare with the first one

      int dist = TMath::Sign(int(pDig1->mX - pDig2->mX), 1) + TMath::Sign(int(pDig1->mY - pDig2->mY), 1); // distance between pads

      if (dist == 1) { // means dig2 is a neighbour of dig1

        if (pDig2->mQ >= pDig1->mQ) {

          iCnt++; // count number of pads with Q more then Q of current pad
        }
      }

    } // second digits loop

    if (iCnt == 0 && mNlocMax < kMaxLocMax) { // this pad has Q more then any neighbour so it's local maximum

      float xStart = o2::hmpid::Digit::lorsX(pDig1->getPadID());

      float yStart = o2::hmpid::Digit::lorsY(pDig1->getPadID());

      float xMin = xStart - param->sizePadX();

      float xMax = xStart + param->sizePadX();

      float yMin = yStart - param->sizePadY();

      float yMax = yStart + param->sizePadY();

      ierflg = fitter->SetParameter(3 * mNlocMax, Form("x%i", mNlocMax), xStart, 0.1, xMin, xMax); // X,Y,Q initial values of the loc max pad

      ierflg = fitter->SetParameter(3 * mNlocMax + 1, Form("y%i", mNlocMax), yStart, 0.1, yMin, yMax); // X, Y constrained to be near the loc max

      ierflg = fitter->SetParameter(3 * mNlocMax + 2, Form("q%i", mNlocMax), pDig1->mQ, 0.1, 0, 10000); // Q constrained to be positive

      mNlocMax++;

    } // if this pad is local maximum

  } // first digits loop

  // Phase 2. Fit loc max number of Mathiesons or add this current cluster to the list

  // case 1 -> no loc max found

  if (mNlocMax == 0) { // case of no local maxima found: pads with same charge...

    mNlocMax = 1;

    mSt = kNoLoc;

    // setClusterParams(mXX, mYY, mCh); //need to fill the AliCluster3D part

    pCluLst->push_back(o2::hmpid::Cluster(*this)); // add new unfolded cluster pCluLst->push_back(o2::hmpid::Cluster(*this));

    pCluLst->back().cleanPointers();

    return mNlocMax;
  }

  // case 2 -> loc max found. Check # of loc maxima

  if (mNlocMax >= kMaxLocMax) {

    LOGP(info, "mNlocMax {} number of digits for cluster {}", mNlocMax, mDigs->size());

    // setClusterParams(mXX, mYY, mCh); // if # of local maxima exceeds kMaxLocMax...

    mSt = kMax;

    pCluLst->push_back(o2::hmpid::Cluster(*this)); //...add this raw cluster

    pCluLst->back().cleanPointers();
  } else { // or resonable number of local maxima to fit and user requested it

    // Now ready for minimization step

    arglist[0] = 500; // number of steps and sigma on pads charges

    arglist[1] = 1.; //

    /*

    ierflg = fitter->ExecuteCommand("SIMPLEX", arglist, 2); // start fitting with Simplex

    if (!ierflg) {

      fitter->ExecuteCommand("MIGRAD", arglist, 2); // fitting improved by Migrad

    }

    if (ierflg) {

      double strategy = 2.;

      ierflg = fitter->ExecuteCommand("SET STR", &strategy, 1); // change level of strategy

      if (!ierflg) {

        ierflg = fitter->ExecuteCommand("SIMPLEX", arglist, 2); // start fitting with Simplex

        if (!ierflg) {

          fitter->ExecuteCommand("MIGRAD", arglist, 2); // fitting improved by Migrad

        }

      }

    }

    */

    double strategy = 2.;

    ierflg = fitter->ExecuteCommand("SET STR", &strategy, 1); // change level of strategy

    if (!ierflg) {

      fitter->ExecuteCommand("MIGRAD", arglist, 2); // fitting improved by Migrad
    }

    if (ierflg) {

      mSt = kAbn; // no convergence of the fit...
    }

    double dummy;

    char sName[80]; // vars to get results from Minuit

    double edm;

    double errdef;

    int nvpar;

    int nparx;

    for (int i = 0; i < mNlocMax; i++) { // store the local maxima parameters

      fitter->GetParameter(3 * i, sName, mXX, mErrX, dummy, dummy); // X

      fitter->GetParameter(3 * i + 1, sName, mYY, mErrY, dummy, dummy); // Y

      fitter->GetParameter(3 * i + 2, sName, mQ, mErrQ, dummy, dummy); // Q

      fitter->GetStats(mChi2, edm, errdef, nvpar, nparx); // get fit infos

      // Printf("********************loc. max. = %i, X= %f, Y = %f, Q = %f**************************",i,mXX,mYY,mQ);

      if (mNlocMax > 1) {

        findClusterSize(i, pSigmaCut); // find clustersize for deconvoluted clusters

        // after this call, fSi temporarly is the calculated size. Later is set again

        // to its original value
      }

      if (mSt != kAbn) {

        if (mNlocMax != 1) {

          mSt = kUnf; // if unfolded
        }

        if (mNlocMax == 1 && mSt != kNoLoc) {

          mSt = kLo1; // if only 1 loc max
        }

        if (!isInPc()) {

          mSt = kEdg; // if Out of Pc
        }

        if (mSt == kNoLoc) {

          mNlocMax = 0; // if with no loc max (pads with same charge..)
        }
      }

      // setClusterParams(mXX, mYY, mCh); //need to fill the AliCluster3D part

      // Printf("********************loc. max. = %i, X= %f, Y = %f, Q = %f**************************",i,mXX,mYY,mQ);

      pCluLst->push_back(o2::hmpid::Cluster(*this)); // add new unfolded cluster

      pCluLst->back().cleanPointers();

      if (mNlocMax > 1) {

        setSize(rawSize); // Original raw size is set again to its proper value
      }
    }
  }

  return mNlocMax;

} // Solve()

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int Cluster::solveMC(std::vector<o2::hmpid::Cluster>* pCluLst, float* pSigmaCut, bool isTryUnfold, std::map<int, std::vector<int>>& resolvedIndicesMap)

{

  // This methode is invoked when the cluster is formed to solve it. Solve the cluster means to try to unfold the cluster

  // into the local maxima number of clusters. This methode is invoked by AliHMPIDRconstructor::Dig2Clu() on cluster by cluster basis.

  // At this point, cluster contains a list of digits, cluster charge and size is precalculated in AddDigit(), position is preset to (-1,-1) in ctor,

  // status is preset to kFormed in AddDigit(), chamber-sector info is preseted to actual values in AddDigit()

  // Method first finds number of local maxima and if it's more then one tries to unfold this cluster into local maxima number of clusters

  // Arguments: pCluLst     - cluster list pointer where to add new cluster(s)

  //            isTryUnfold - flag to switch on/off unfolding

  //   Returns: number of local maxima of original cluster

  const auto param = o2::hmpid::Param::instanceNoGeo();

  if (!mDigs) {

    LOGP(fatal, "digits are missing in the cluster");
  }

  const int kMaxLocMax = 6; // max allowed number of loc max for fitting

  coG(); // First calculate CoG for the given cluster

  int iCluCnt = pCluLst->size(); // get current number of clusters already stored in the list by previous operations

  int rawSize = mSi; // get current raw cluster size

  if (rawSize > 100) {

    mSt = kBig;
  } else if (isTryUnfold == false) {

    mSt = kNot;
  } else if (rawSize == 1) {

    mSt = kSi1;
  }

  if (rawSize > 100 || isTryUnfold == false || rawSize == 1) { // No deconv if: 1 - big cluster (also avoid no zero suppression!)

    // setClusterParams(mXX, mYY, mCh); //                               2 - flag is set to FALSE

    // new ((*pCluLst)[iCluCnt++]) Cluster(*this); //                      3 - size = 1

    pCluLst->push_back(o2::hmpid::Cluster(*this));

    LOGP(info, "mNlocMax :: no deconv number of digits for cluster {}", mDigs->size());

    // if(useMC) { // EF: TODO: Set useMC cond here ?

    pCluLst->back().setInfoFromDigits(pCluLst->size());

    // pCluLst->back().setMC();

    // , pCluLst->size() to set the index of the cluster in teh vector to each element

    // if(useMC) {

    //  pCluLst->back().setDigitTruth();

    //}

    pCluLst->back().cleanPointers();

    return 1; // add this raw cluster
  }

  //  Phase 0. Initialise Fitter

  double arglist[10]{0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  float ierflg = 0.;

  TVirtualFitter* fitter = TVirtualFitter::Fitter((TObject*)this, 3 * 6); // initialize Fitter

  if (fitter == nullptr) {

    LOG(fatal) << "TVirtualFitter could not be created";

    return 1;
  }

  arglist[0] = -1;

  ierflg = fitter->ExecuteCommand("SET PRI", arglist, 1); // no printout

  ierflg = fitter->ExecuteCommand("SET NOW", arglist, 0); // no warning messages

  arglist[0] = 1;

  ierflg = fitter->ExecuteCommand("SET GRA", arglist, 1); // force Fitter to use my gradient

  fitter->SetFCN(Cluster::fitFunc);

  // Phase 1. Find number of local maxima. Strategy is to check if the current pad has QDC more then all neigbours. Also find the box contaning the cluster

  mNlocMax = 0;

  for (int iDig1 = 0; iDig1 < rawSize; iDig1++) { // first digits loop

    auto pDig1 = (*mDigs)[iDig1]; // take next digit

    int iCnt = 0; // counts how many neighbouring pads has QDC more then current one

    for (int iDig2 = 0; iDig2 < rawSize; iDig2++) { // loop on all digits again

      if (iDig1 == iDig2) {

        continue;

      } // the same digit, no need to compare

      auto pDig2 = (*mDigs)[iDig2]; // take second digit to compare with the first one

      int dist = TMath::Sign(int(pDig1->mX - pDig2->mX), 1) + TMath::Sign(int(pDig1->mY - pDig2->mY), 1); // distance between pads

      if (dist == 1) { // means dig2 is a neighbour of dig1

        if (pDig2->mQ >= pDig1->mQ) {

          iCnt++; // count number of pads with Q more then Q of current pad
        }
      }

    } // second digits loop

    if (iCnt == 0 && mNlocMax < kMaxLocMax) { // this pad has Q more then any neighbour so it's local maximum

      float xStart = o2::hmpid::Digit::lorsX(pDig1->getPadID());

      float yStart = o2::hmpid::Digit::lorsY(pDig1->getPadID());

      float xMin = xStart - param->sizePadX();

      float xMax = xStart + param->sizePadX();

      float yMin = yStart - param->sizePadY();

      float yMax = yStart + param->sizePadY();

      ierflg = fitter->SetParameter(3 * mNlocMax, Form("x%i", mNlocMax), xStart, 0.1, xMin, xMax); // X,Y,Q initial values of the loc max pad

      ierflg = fitter->SetParameter(3 * mNlocMax + 1, Form("y%i", mNlocMax), yStart, 0.1, yMin, yMax); // X, Y constrained to be near the loc max

      ierflg = fitter->SetParameter(3 * mNlocMax + 2, Form("q%i", mNlocMax), pDig1->mQ, 0.1, 0, 10000); // Q constrained to be positive

      mNlocMax++;

    } // if this pad is local maximum

  } // first digits loop

  // Phase 2. Fit loc max number of Mathiesons or add this current cluster to the list

  // case 1 -> no loc max found

  if (mNlocMax == 0) { // case of no local maxima found: pads with same charge...

    LOGP(info, "mNlocMax {} number of digits for cluster {}", mNlocMax, mDigs->size());

    mNlocMax = 1;

    mSt = kNoLoc;

    // setClusterParams(mXX, mYY, mCh); //need to fill the AliCluster3D part

    pCluLst->push_back(o2::hmpid::Cluster(*this)); // add new unfolded cluster pCluLst->push_back(o2::hmpid::Cluster(*this));

    pCluLst->back().setInfoFromDigits(pCluLst->size());

    // pCluLst->back().setMC();

    pCluLst->back().cleanPointers();

    return mNlocMax;
  }

  // case 2 -> loc max found. Check # of loc maxima

  if (mNlocMax >= kMaxLocMax) {

    LOGP(info, "mNlocMax {} number of digits for cluster {}", mNlocMax, mDigs->size());

    // setClusterParams(mXX, mYY, mCh); // if # of local maxima exceeds kMaxLocMax...

    mSt = kMax;

    pCluLst->push_back(o2::hmpid::Cluster(*this)); //...add this raw cluster

    pCluLst->back().setInfoFromDigits(pCluLst->size());

    // pCluLst->back().setMC();

    pCluLst->back().cleanPointers();
  } else { // or resonable number of local maxima to fit and user requested it

    // Now ready for minimization step

    arglist[0] = 500; // number of steps and sigma on pads charges

    arglist[1] = 1.; //

    /*

    ierflg = fitter->ExecuteCommand("SIMPLEX", arglist, 2); // start fitting with Simplex

    if (!ierflg) {

      fitter->ExecuteCommand("MIGRAD", arglist, 2); // fitting improved by Migrad

    }

    if (ierflg) {

      double strategy = 2.;

      ierflg = fitter->ExecuteCommand("SET STR", &strategy, 1); // change level of strategy

      if (!ierflg) {

        ierflg = fitter->ExecuteCommand("SIMPLEX", arglist, 2); // start fitting with Simplex

        if (!ierflg) {

          fitter->ExecuteCommand("MIGRAD", arglist, 2); // fitting improved by Migrad

        }

      }

    }

    */

    double strategy = 2.;

    ierflg = fitter->ExecuteCommand("SET STR", &strategy, 1); // change level of strategy

    if (!ierflg) {

      fitter->ExecuteCommand("MIGRAD", arglist, 2); // fitting improved by Migrad
    }

    if (ierflg) {

      mSt = kAbn; // no convergence of the fit...
    }

    double dummy;

    char sName[80]; // vars to get results from Minuit

    double edm;

    double errdef;

    int nvpar;

    int nparx;

    for (int i = 0; i < mNlocMax; i++) { // store the local maxima parameters

      fitter->GetParameter(3 * i, sName, mXX, mErrX, dummy, dummy); // X

      fitter->GetParameter(3 * i + 1, sName, mYY, mErrY, dummy, dummy); // Y

      fitter->GetParameter(3 * i + 2, sName, mQ, mErrQ, dummy, dummy); // Q

      fitter->GetStats(mChi2, edm, errdef, nvpar, nparx); // get fit infos

      // Printf("********************loc. max. = %i, X= %f, Y = %f, Q = %f**************************",i,mXX,mYY,mQ);

      if (mNlocMax > 1) {

        std::vector<int> indicesResolved;
        findClusterSizeMC(i, pSigmaCut, indicesResolved); // find clustersize for deconvoluted clusters

        // add "local" indices to map of clusterNumer i
        resolvedIndicesMap[i] = indicesResolved;

        findClusterSize(i, pSigmaCut); // find clustersize for deconvoluted clusters

        // after this call, fSi temporarly is the calculated size. Later is set again

        // to its original value

        /*if(useMC) {

          findDigitsForCluster(i, solve);

        }*/
      }

      if (mSt != kAbn) {

        if (mNlocMax != 1) {

          mSt = kUnf; // if unfolded
        }

        if (mNlocMax == 1 && mSt != kNoLoc) {

          mSt = kLo1; // if only 1 loc max
        }

        if (!isInPc()) {

          mSt = kEdg; // if Out of Pc
        }

        if (mSt == kNoLoc) {

          mNlocMax = 0; // if with no loc max (pads with same charge..)
        }
      }

      // setClusterParams(mXX, mYY, mCh); //need to fill the AliCluster3D part

      // Printf("********************loc. max. = %i, X= %f, Y = %f, Q = %f**************************",i,mXX,mYY,mQ);

      LOGP(info, "cluNumber {} mNlocMax {} number of digits for cluster {}", i, mNlocMax, mDigs->size());

      pCluLst->push_back(o2::hmpid::Cluster(*this)); // add new unfolded cluster

      pCluLst->back().setInfoFromDigits(pCluLst->size());

      // pCluLst->back().setMC();

      pCluLst->back().cleanPointers();

      if (mNlocMax > 1) {

        setSize(rawSize); // Original raw size is set again to its proper value
      }
    }
  }

  return mNlocMax;

} // Solve()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Estimate of the clustersize for a deconvoluted cluster

void Cluster::findClusterSizeMC(int i, float* pSigmaCut, std::vector<int>& indicesResolved)

{

  // std::vector<int> indexResolved;

  // auto indexUnresolved = getUnresolvedIndexes();

  int size = 0;

  for (int iDig = 0; iDig < mSi; iDig++) { // digits loop

    auto pDig = dig(iDig); // take digit

    int iCh = pDig->mCh;

    double qPad = mQ * o2::hmpid::Digit::intMathieson(x(), y(), pDig->getPadID()); // pad charge  pDig->

    //  AliDebug(1,Form("Chamber %i X %i Y %i SigmaCut %i pad %i qpadMath %8.2f qPadRaw %8.2f Qtotal %8.2f cluster n.%i",

    //                 iCh, o2::hmpid::Digit::a2X(pDig->getPadID()), o2::hmpid::Digit::a2Y(pDig->getPadID()),

    //                 pSigmaCut[iCh],iDig,qPad,pDig->mQ,mQRaw,i));

    if (qPad > pSigmaCut[iCh]) {

      indicesResolved.push_back(iDig); // ef : added to track indexes of resolved clusters
      LOGP(info, "findClusterSizeMC :: added local iDig {} ", iDig);

      size++;
    }
  }

  LOGP(info, "findClusterSizeMC {}", size);

  //  AliDebug(1,Form(" Calculated size %i",size));

  if (size > 0) {

    setSize(size); // in case of size == 0, original raw clustersize used
  }

  // ef : set the resolved indexes

  // done above instead

  // this->setResolvedIndexes(indexResolved);

  // setResolvedIndex(indexResolved);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Estimate of the clustersize for a deconvoluted cluster

void Cluster::findClusterSize(int i, float* pSigmaCut)

{

  int size = 0;

  for (int iDig = 0; iDig < mSi; iDig++) { // digits loop

    auto pDig = dig(iDig); // take digit

    int iCh = pDig->mCh;

    double qPad = mQ * o2::hmpid::Digit::intMathieson(x(), y(), pDig->getPadID()); // pad charge  pDig->

    //  AliDebug(1,Form("Chamber %i X %i Y %i SigmaCut %i pad %i qpadMath %8.2f qPadRaw %8.2f Qtotal %8.2f cluster n.%i",

    //                 iCh, o2::hmpid::Digit::a2X(pDig->getPadID()), o2::hmpid::Digit::a2Y(pDig->getPadID()),

    //                 pSigmaCut[iCh],iDig,qPad,pDig->mQ,mQRaw,i));

    if (qPad > pSigmaCut[iCh]) {
      size++;
    }
  }

  //  AliDebug(1,Form(" Calculated size %i",size));

  if (size > 0) {
    setSize(size); // in case of size == 0, original raw clustersize used
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t Cluster::isInPc()

{

  // Check if (X,Y) position is inside the PC limits

  // Arguments:

  //   Returns: True or False

  const auto param = o2::hmpid::Param::instanceNoGeo();

  if (!mDigs) {

    LOGP(fatal, "digits are missing in the cluster");
  }

  int pc = (*mDigs)[0]->getPh(); // (o2::hmpid::Digit*)&mDigs.at(iDig)

  if (mXX < param->minPcX(pc) || mXX > param->maxPcX(pc) || mYY < param->minPcY(pc) || mYY > param->maxPcY(pc)) {

    return false;
  }

  return true;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Cluster::digAdd(const Digit* pDig)

{

  // Adds a given digit to the list of digits belonging to this cluster, cluster is not owner of digits

  // Arguments: pDig - pointer to digit to be added

  // Returns: none

  if (!mDigs) {

    LOGP(fatal, "digits are not set to the cluster");
  }

  if (mDigs->size() == 0) { // create list of digits in the first invocation

    mSi = 0;
  }

  // fDigs->Add(pDig);

  mDigs->push_back(pDig);

  mSt = kFrm;

  mSi++;

  return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Cluster::reset()

{

  //

  cleanPointers();

  // mDigs={0x0};

  mSt = kEmp;

  mQRaw = mQ = 0;

  mXX = mYY = 0;

  mCh = mSi = mBox = mNlocMax = mMaxQpad = -1;

  mMaxQ = mErrQ = mErrX = mErrY = mChi2 = -1; // empty ctor
}

} // namespace hmpid

} // namespace o2
