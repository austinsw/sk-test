#ifndef SK_H
#define SK_H

//--------------------------------------------------------------------------------------------------//
//  DFTB+: general package for performing fast atomistic simulations                                //
//  Copyright (C) 2006 - 2022  DFTB+ developers group                                               //
//                                                                                                  //
//  See the LICENSE file for terms of usage and distribution.                                       //
//--------------------------------------------------------------------------------------------------//

#include <iostream>
#include "t_orbitals.h"
#include <assert.h>
#include <vector>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

// Contains code to perform the sk rotations of matrix elements from the parameterization
// orientation along <0,0,1> to the one needed in the actual calculation.
//
// To do: Transformations to give the derivatives with respect to ll, mm and nn.
// Base on "Compact expression for the angular dependence of tight-binding hamiltonian matrix
// elements", A. V. Podolskiy and P. Vogl, Phys. Rev.  B 69, 233101 (2004).
//
// Caveat: Only angular momenta up to f are currently allowed
class dftbp_dftb_sk {

  // Maximal angular momentum, for which rotations are present
  int mAngRot_ = 3;

  // A temporary replacement for Fortran's built in maxval function
  int maxval(vector<vector<int>> const &vv) {
    vector<int> tmp;
    for (int i = 0; i < vv.size(); i++)
      tmp.insert(tmp.end(),vv[i].begin(),vv[i].end());
    return *max_element(tmp.begin(),tmp.end());
  }

  // Driver for making the non-SCC hhamiltonian or overlap matrices for a given diatomic block
  // Caveat: Only angular momenta up to f are currently allowed
  void rotateH0(

    // the rectangular matrix containing the resulting diatomic matrix elements
    vector<vector<double>> &hh,

    // Slater-Koster table for dimer of species i-j
    vector<double> const &skIntegs,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Chemical species of atom i
    int const &iSp1,

    // chemical species of atom j
    int const &iSp2,
    
    // Information about the orbitals of chemical species in the system.
    TOrbitals const &orb

  ) 
  {
    int iCol, iRow, ind, iSh1, iSh2;
    int ang1, ang2, nOrb1, nOrb2;
    vector<double> *pSK;
    //double tmpH[2*mAngRot_+1][2*mAngRot_+1];
    vector<vector<double>> tmpH(2*mAngRot_+1, vector<double> (2*mAngRot_+1, 0.0));
    double sign;

    assert(maxval(orb.angShell) <=  mAngRot_);
    assert(hh.size() >= orb.nOrbSpecies[iSp1] && hh[0].size() >= orb.nOrbSpecies[iSp2]);

    fill(hh.begin(), hh.end(), 0.0); //hh = 0.0;
    ind = 1;
    iCol = 1;
    for (iSh1 = 1; iSh1 <= orb.nShell[iSp1]; iSh1++) {
      ang1 = orb.angShell[iSh1][iSp1];
      nOrb1 = 2 * ang1 + 1;
      iRow = 1;
      for (iSh2 = 1; iSh2 <= orb.nShell[iSp2]; iSh2++) {
        ang2 = orb.angShell[iSh2][iSp2];
        nOrb2 = 2 * ang2 + 1;
        assert(skIntegs.size() >= ind + min(ang1, ang2));
        pSK = &vector<double> (skIntegs.begin() + ind, skIntegs.begin() + ind + min(ang1, ang2));
        switch (ang1) {
          case 0:
            switch (ang2) {
              case 0:
                ss(tmpH,*pSK);
              case 1:
                sp(tmpH,ll,mm,nn,*pSK);
              case 2:
                sd(tmpH,ll,mm,nn,*pSK);
              case 3:
                sf(tmpH,ll,mm,nn,*pSK);
            }
          case 1:
            switch (ang2) {
              case 0:
                sp(tmpH,ll,mm,nn,*pSK);
              case 1:
                pp(tmpH,ll,mm,nn,*pSK);
              case 2:
                pd(tmpH,ll,mm,nn,*pSK);
              case 3:
                pf(tmpH,ll,mm,nn,*pSK);
            }
          case 2:
            switch (ang2) {
              case 0:
                sd(tmpH,ll,mm,nn,*pSK);
              case 1:
                pd(tmpH,ll,mm,nn,*pSK);
              case 2:
                dd(tmpH,ll,mm,nn,*pSK);
              case 3:
                df(tmpH,ll,mm,nn,*pSK);
            }
          case 3:
            switch (ang2) {
              case 0:
                sf(tmpH,ll,mm,nn,*pSK);
              case 1:
                pf(tmpH,ll,mm,nn,*pSK);
              case 2:
                df(tmpH,ll,mm,nn,*pSK);
              case 3:
                ff(tmpH,ll,mm,nn,*pSK);
            }
        }

        if (ang1 <= ang2) {
          //hh[iRow:iRow+nOrb2-1][iCol:iCol+nOrb1-1] = tmpH[1:nOrb2][1:nOrb1];
          for (int i = 0; i < nOrb2; i++)
            for (int j = 0; j < nOrb1; j++)
              hh[i+iRow-1][j+iCol-1] = tmpH[i][j];
        } else {
          //hh[iRow:iRow+nOrb2-1][iCol:iCol+nOrb1-1] = (-1.0)**(ang1+ang2) * transpose(tmpH(1:nOrb1,1:nOrb2))
          sign = pow(-1.0, (ang1+ang2));
          for (int i = 0; i < nOrb2; i++)
            for (int j = 0; j < nOrb1; j++)
              hh[i+iRow-1][j+iCol-1] = sign * tmpH[i][j];
        }
        ind = ind + min(ang1,ang2) + 1;
        iRow = iRow + nOrb2;
      }
      iCol = iCol + nOrb1;
    }
  }

  // rotation routine for interaction of an s orbital with an s orbital
  void ss(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 1);
    assert(hh.size() >= 1 && hh[0].size() >= 1);
    hh[0][0] = sk[0];
  }

  // rotation routine for interaction of an s orbital with a p orbital
  void sp(
    
    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 1);
    assert(hh.size() >= 3 && hh[0].size() >= 1);

    hh[1][1] = mm * sk[1];
    hh[2][1] = nn * sk[1];
    hh[3][1] = ll * sk[1];
  }

  // rotation routine for interaction of an s orbital with a d orbital
  void sd(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 1);
    assert(hh.size() >= 5 && hh[0].size() >= 1);

    hh[1][1] = ll * mm * sqrt(3.0) * sk[1];
    hh[2][1] = mm * sqrt(3.0) * nn * sk[1];
    hh[3][1] = (3.0 / 2.0 * pow(nn, 2) - 1.0 / 2.0) * sk[1];
    hh[4][1] = ll * sqrt(3.0) * nn * sk[1];
    hh[5][1] = (2.0 * ll * 2 - 1.0 + pow(nn, 2)) * sqrt(3.0) * sk[1] / 2.0;
  }

  // rotation routine for interaction of an s orbital with an f orbital
  void sf(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 1);
    assert(hh.size() >= 1 && hh[0].size() >= 7);

    hh[1][1] = sqrt(2.0) * mm * (4.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(5.0) * sk[1] / 4.0;
    hh[2][1] = ll * mm * sqrt(15.0) * nn * sk[1];
    hh[3][1] = sqrt(2.0) * mm * sqrt(3.0) * (5.0 * pow(nn, 2) - 1.0) * sk[1] / 4.0;
    hh[4][1] = (nn * (5.0 * pow(nn, 2) - 3.0) * sk[1]) / 2.0;
    hh[5][1] = sqrt(2.0) * ll * sqrt(3.0) * (5.0 * pow(mm, 2) - 1.0) * sk[1] / 4.0;
    hh[6][1] = (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(15.0) * nn * sk[1] / 2.0;
    hh[7][1] = sqrt(2.0) * ll * (4.0 * pow(ll, 2) - 3.0 + 3.0 * pow(nn, 2)) * sqrt(5.0) * sk[1] / 4.0;
  }

  // rotation routine for interaction of a p orbital with a p orbital
  void pp(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 2);
    assert(hh.size() >= 3 && hh[0].size() >= 3);

    hh[1][1] = (1.0 - pow(nn, 2) - pow(ll, 2)) * sk[1] + (pow(nn, 2) + pow(ll, 2)) * sk[2];
    hh[2][1] = nn * mm * sk[1] - nn * mm * sk[2];
    hh[3][1] = ll * mm * sk[1] - ll * mm * sk[2];
    hh[1][2] = hh[2][1];
    hh[2][2] = pow(nn, 2) * sk[1] + (1.0 - pow(nn, 2)) * sk[2];
    hh[3][2] = nn * ll * sk[1] - nn * ll * sk[2];
    hh[1][3] = hh[3][1];
    hh[2][3] = hh[3][2];
    hh[3][3] = pow(ll, 2) * sk[1] + (1.0 - pow(ll, 2)) * sk[2];
  }

  // rotation routine for interaction of a p orbital with a d orbital
  void pd(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 2);
    assert(hh.size() >= 3 && hh[0].size() >= 5);

    hh[1][1] = -(-1.0 + pow(nn, 2) + pow(ll, 2)) * ll * sqrt(3.0) * sk[1] + ((2.0 * pow(nn,2) + 2.0 * pow(ll, 2) -1.0) * ll * sk[2]);
    hh[2][1] = -(-1.0 + pow(nn, 2) + pow(ll, 2)) * sqrt(3.0) * nn * sk[1] + ((2.0 * pow(nn,2) + 2.0 * pow(ll, 2) -1.0) * nn * sk[2]);
    hh[3][1] = mm * (3.0 * pow(nn, 2) -1.0) * sk[1] / 2.0 - sqrt(3.0) * pow(nn, 2) * mm * sk[2];
    hh[4][1] = mm * ll * sqrt(3.0) * nn * sk[1] - 2.0 * ll * mm * nn * sk[2];
    hh[5][1] = mm * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(3.0) * sk[1] / 2.0 - (pow(nn, 2) + 2.0) * pow(ll, 2) * mm * sk[2];
    hh[1][2] = ll * mm * nn * sqrt(3.0) * sk[1] - 2.0 * nn * ll * mm * sk[2];
    hh[2][2] = mm * pow(nn, 2) * sqrt(3.0) * sk[1] - (2.0 * pow(nn, 2) - 1.0) * mm * sk[2];
    hh[3][2] = (nn * (3.0 * pow(nn, 2) -1.0) * sk[1]) / 2.0 - nn * sqrt(3.0) * (-1.0 + pow(nn, 2)) * sk[2];
    hh[4][2] = ll * pow(nn, 2) * sqrt(3.0) * sk[1] - (2.0 * pow(nn, 2) - 1.0) * ll * sk[2];
    hh[5][2] = (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * nn * sqrt(3.0) * sk[1] / 2.0 - (nn * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sk[2]);
    hh[1][3] = pow(ll, 2) * mm * sqrt(3.0) * sk[1] - (2.0 * pow(ll, 2) - 1.0) * mm * sk[2];
    hh[2][3] = ll * mm * sqrt(3.0) * nn * sk[1] - 2.0 * mm * ll * nn * sk[2];
    hh[3][3] = (ll * (3.0 * pow(nn, 2) - 1.0) * sk[1]) / 2.0 - sqrt(3.0) * pow(nn, 2) * ll * sk[2];
    hh[4][3] = pow(ll, 2) * sqrt(3.0) * nn * sk[1] - (2.0 * pow(ll, 2) - 1.0) * nn * sk[2];
    hh[5][3] = ll * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(3.0) * sk[1] / 2.0 - ((pow(nn, 2) -2.0 + 2.0 * pow(ll, 2)) * ll * sk[2]);
  }

  // rotation routine for interaction of a p orbital with an f orbital
  void pf(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 2);
    assert(hh.size() >= 3 && hh[0].size() >= 7);

    hh[1][1] = -(-1.0 + pow(nn, 2) + pow(ll, 2)) * (4.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(2.0) * sqrt(5.0) * sk[1] / 4.0 
              + sqrt(15.0) * (pow(nn, 4) - pow(nn, 2) + 5.0 * pow(nn, 2) * pow(ll, 2) - 3.0 * pow(ll, 2) + 4.0 * pow(ll, 4)) * sk[2] / 4.0;
    hh[2][1] = -(-1.0 + pow(nn, 2) + pow(ll, 2)) * ll * sqrt(15.0) * nn * sk[1] + (3.0 * pow(nn, 2) + 3.0 * pow(ll, 2) - 2.0) * ll * nn * sqrt(10.0) * sk[2] / 2.0;
    hh[3][1] = -(-1.0 + pow(nn, 2) + pow(ll, 2)) * sqrt(2.0) * sqrt(3.0) * (5.0 * pow(nn, 2) - 1.0) * sk[1] / 4.0 
              + (15.0 / 4.0 * pow(nn, 4) + 15.0 / 4.0 * pow(ll, 2) * pow(nn, 2) - 11.0 / 4.0 * pow(nn, 2) - pow(ll, 2) / 4.0) * sk[2];
    hh[4][1] = mm * nn * (5.0 * pow(nn, 2) - 3.0) * sk[1] / 2.0 - (5.0 * pow(nn, 2) -1.0) * sqrt(3.0) * sqrt(2.0) * nn * mm * sk[2] / 4.0;
    hh[5][1] = mm * ll * sqrt(2.0) * sqrt(3.0) * (5.0 * pow(nn, 2) - 1.0) * sk[1] / 4.0 - (15.0 * pow(nn, 2) - 1.0) * ll * mm * sk[2] / 4.0;
    hh[6][1] = mm * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(15.0) * nn * sk[1] / 2.0 - (3.0 * pow(nn, 2) + 6.0 * pow(ll, 2) - 1.0) * nn * mm * sqrt(10.0) * sk[2] / 4.0;
    hh[7][1] = mm * ll * (4.0) * pow(ll, 2) - 3.0 + 3.0 * pow(nn, 2) * sqrt(2.0) * sqrt(5.0) * sk[1] / 4.0 - ll * mm * sqrt(15.0)
              * (3.0 * pow(nn, 2) + 4.0 * pow(ll, 2) -1.0) * sk[2] / 4.0;
    hh[1][2] = sqrt(2.0) * mm * (4.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * nn * sqrt(5.0) * sk[1] / 4.0 - mm * (4.0 * pow(ll, 2) - 1.0 + pow(nn, 2))
              * sqrt(15.0) * nn * sk[2] / 4.0;
    hh[2][2] = ll * mm * pow(nn, 2) * sqrt(15.0) * sk[1] - (3.0 * pow(nn, 2) - 1.0) * ll * mm * sqrt(10.0) * sk[2] / 2.0;
    hh[3][2] = sqrt(2.0) * mm * nn * sqrt(3.0) * (5.0 * pow(nn, 2) - 1.0) * sk[1] / 4.0 - (15.0 * pow(nn, 2) -11.0) * nn * mm * sk[2] / 4.0;
    hh[4][2] = (pow(nn, 2) * (5.0 * pow(nn, 2) -3.0) * sk[1]) / 2.0 - (5.0 * pow(nn, 2) - 1.0) * sqrt(3.0) * sqrt(2.0) 
              * (-1.0 + pow(nn, 2)) * sk[2] / 4.0;
    hh[5][2] = sqrt(2.0) * ll * nn * sqrt(3.0) * (5.0 * pow(nn, 2) - 1.0) * sk[1] / 4.0 - (15.0 * pow(nn, 2) - 11.0) * nn * ll * sk[2] / 4.0;
    hh[6][2] = (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * pow(nn, 2) * sqrt(15.0) * sk[1] / 2.0 - (3.0 * pow(nn, 2) - 1.0)
              * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(10.0) * sk[2] / 4.0;
    hh[7][2] = sqrt(2.0) * ll * (4.0 * pow(ll, 2) - 3.0 + 3.0 * pow(nn, 2)) * nn * sqrt(5.0) * sk[1] / 4.0
              - ll * (4.0 * pow(ll, 2) - 3.0 + 3.0 * pow(nn, 2)) * sqrt(15.0) * nn * sk[2] / 4.0;
    hh[1][3] = ll * mm * (4.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(2.0) * sqrt(5.0) * sk[1] / 4.0 - ll * mm * sqrt(15.0)
              * (pow(nn, 2) + 4.0 * pow(ll, 2) - 3.0) * sk[2] / 4.0;
    hh[2][3] = pow(ll, 2) * mm * sqrt(15.0) * nn * sk[1] - (3.0 * pow(ll, 2) -1.0) * nn * mm * sqrt(10.0) * sk[2] / 2.0;
    hh[3][3] = ll * mm * sqrt(2.0) * sqrt(3.0) * (5.0 * pow(nn, 2) - 1.0) * sk[1] / 4.0 - (15.0 * pow(nn, 2) - 1.0) * ll * mm * sk[2] / 4.0;
    hh[4][3] = (ll * nn * (5.0 * pow(nn, 2) -3.0) * sk[1]) / 2.0 - (5.0 * pow(nn, 2) - 1.0) * sqrt(3.0) * sqrt(2.0) * nn * ll * sk[2] / 4.0;
    hh[5][3] = pow(ll, 2) * sqrt(2.0) * sqrt(3.0) * (5.0 * pow(nn, 2) - 1.0) * sk[1] / 4.0 + (-15.0 / 4.0 * pow(ll, 2) * pow(nn, 2) 
              + 5.0 / 4.0 * pow(nn, 2) - 1.0 / 4.0 + pow(ll, 2) / 4.0) * sk[2];
    hh[6][3] = ll * (2.0 * pow(ll, 2) -1.0 + pow(nn, 2)) * sqrt(15.0) * nn * sk[1] / 2.0 - (3.0 * pow(nn, 2) + 6.0 * pow(ll, 2) - 5.0)
              * nn * ll * sqrt(10.0) * sk[2] / 4.0;
    hh[7][3] = pow(ll, 2) * (4.0 * pow(ll, 2) - 3.0 + 3.0 * pow(nn, 2)) * sqrt(2.0) * sqrt(5.0) * sk[1] / 4.0 - sqrt(15.0) 
              * (3.0 * pow(ll, 2) * pow(nn, 2) - pow(nn, 2) - 5.0 * pow(ll, 2) + 4.0 * pow(ll, 4) + 1.0) * sk[2] / 4.0;

  }

  // rotation routine for interaction of a d orbital with a d orbital
  void dd(
    
    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 3);
    assert(hh.size() >= 5 && hh[0].size() >= 5);

    hh[1][1] = -3.0 * pow(ll, 2) * (-1.0 + pow(nn, 2) + pow(ll, 2)) * sk[1] + (4.0 * pow(ll, 2) * pow(nn, 2) - pow(nn, 2) + 4.0 
              * pow(ll, 4) - 4.0 * pow(ll, 2) + 1.0) * sk[2] + (-pow(ll, 2) * pow(nn, 2) + pow(nn, 2) + pow(ll, 2) - pow(ll, 4)) * sk[3];
    hh[2][1] = -3.0 * ll * (-1.0 + pow(nn, 2) + pow(ll, 2)) * nn * sk[1] + (4.0 * pow(nn, 2) + 4.0 * pow(ll, 2) - 3.0) * nn * ll * sk[2]
              - ll * (pow(nn, 2) + pow(ll, 2)) * nn * sk[3];
    hh[3][1] = ll * mm * sqrt(3.0) * (3.0 * pow(nn, 2) - 1.0) * sk[1] / 2.0 - 2.0 * sqrt(3.0) * mm * ll * pow(nn, 2) * sk[2]
              + ll * mm * (pow(nn, 2) + 1.0) * sqrt(3.0) * sk[3] / 2.0;
    hh[4][1] = 3.0 * pow(ll, 2) * mm * nn * sk[1] - (4.0 * pow(ll, 2) - 1.0) * nn * mm * sk[2] + mm * (-1.0 + pow(ll,2)) * nn * sk[3];
    hh[5][1] = 3.0 / 2.0 * mm * ll * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sk[1] - 2.0 * mm * ll * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2))
              * sk[2] + mm * ll * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sk[3] / 2.0;
    hh[1][2] = hh[2][1];
    hh[2][2] = -3.0 * (-1.0 + pow(nn, 2) + pow(ll, 2)) * pow(nn, 2) * sk[1] + (4.0 * pow(nn, 4) - 4.0 * pow(nn, 2) + 4.0 * pow(ll, 2)
              * pow(nn, 2) + 1.0 - pow(ll, 2)) * sk[2] -(-1.0 + nn) * (pow(nn, 3) + pow(nn, 2) + pow(ll, 2)) * sk[3];
    hh[3][2] = mm * sqrt(3.0) * nn * (3.0 * pow(nn, 2) - 1.0) * sk[1] / 2.0 - nn * sqrt(3.0) * (2.0 * pow(nn, 2) -1.0) * mm * sk[2]
              + (-1.0 + pow(nn, 2)) * sqrt(3.0) * nn * mm * sk[3] / 2.0;
    hh[4][2] = 3.0 * mm * ll * pow(nn, 2) * sk[1] - (4.0 * pow(nn, 2) - 1.0) * mm * ll * sk[2] + ll * mm * (-1.0 + pow(nn, 2)) * sk[3];
    hh[5][2] = 3.0 / 2.0 * mm * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * nn * sk[1] - (2.0 * pow(nn, 2) - 1.0 + 4.0 * pow(ll, 2)) 
              * nn * mm * sk[2] + mm * (pow(nn, 2) + 2.0 * pow(ll, 2) + 1.0) * nn * sk[3] / 2.0;
    hh[1][3] = hh[3][1];
    hh[2][3] = hh[3][2];
    hh[3][3] = pow(3.0 * pow(nn, 2) - 1.0, 2) * sk[1] / 4.0 - 3.0 * (-1.0 + pow(nn, 2)) * pow(nn, 2) * sk[2] + 3.0 / 4.0 * pow(-1.0 + pow(nn, 2), 2) * sk[3];
    hh[4][3] = ll * (3.0 * pow(nn, 2) - 1.0) * sqrt(3.0) * nn * sk[1] / 2.0 - (2.0 * pow(nn, 2) - 1.0) * ll * nn * sqrt(3.0) * sk[2]
              + nn * ll * sqrt(3.0) * (-1.0 + pow(nn, 2)) * sk[3] / 2.0;
    hh[5][3] = (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * (3.0 * pow(nn, 2) - 1.0) * sqrt(3.0) * sk[1] / 4.0 - (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2))
              * pow(nn, 2) * sqrt(3.0) * sk[2] + sqrt(3.0) * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * (pow(nn, 2) + 1.0) * sk[3] / 4.0;
    hh[1][4] = hh[4][1];
    hh[2][4] = hh[4][2];
    hh[3][4] = hh[4][3];
    hh[4][4] = 3.0 * pow(ll, 2) * pow(nn, 2) * sk[1] + (-4.0 * pow(ll, 2) * pow(nn, 2) + pow(nn, 2) + pow(ll, 2)) * sk[2] + (-1.0 + nn)
              * (-nn + pow(ll, 2) * nn - 1.0 + pow(ll, 2)) * sk[3];
    hh[5][4] = 3.0 / 2.0 * ll * (2.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * nn * sk[1] - ((2.0 * pow(nn, 2) - 3.0 + 4.0 * pow(ll, 2)) * nn
              * ll * sk[2]) + (ll * (pow(nn, 2) - 3.0 + 2.0 * pow(ll, 2)) * nn * sk[3]) / 2.0;
    hh[1][5] = hh[5][1];
    hh[2][5] = hh[5][2];
    hh[3][5] = hh[5][3];
    hh[4][5] = hh[5][4];
    hh[5][5] = 3.0 / 4.0 * (pow(2.0 * pow(ll, 2) - 1.0 + pow(nn, 2), 2)) * sk[1] + ((-pow(nn, 4) + pow(nn, 2) - 4.0 * pow(ll, 2) * pow(nn, 2)
              - 4.0 * pow(ll, 4) + 4.0 * pow(ll, 2)) * sk[2]) + (pow(nn, 4) / 4.0 + (pow(ll, 2) * pow(nn, 2)) + pow(nn, 2) / 2.0 + 1.0 / 4.0
              - pow(ll, 2) + pow(ll, 4) * sk[3]);
  }

  // rotation routine for interaction of a d orbital with an f orbital
  void df(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 3);
    assert(hh.size() >= 5 && hh[0].size() >= 7);

    hh[1][1] = -ll *(-1.0 + pow(nn, 2) + pow(ll, 2)) * (4.0 * pow(ll, 2) - 1.0 + pow(nn, 2)) * sqrt(6.0) * sqrt(5.0) * sk[1] / 4.0
              + sqrt(15.0) * ll * (2.0 * pow(nn, 4) - 5.0 * pow(nn, 2) + 10.0 * pow(ll, 2) * pow(nn, 2) + 3.0 - 10.0 * pow(ll, 2)
              + 8.0 * pow(ll, 4)) * sk[2] / 4.0 - sqrt(6.0) * ll * (pow(nn, 4) + 5.0 * pow(ll, 2) - 4.0 * pow(nn, 2) + 1.0 + 4.0
              * pow(ll, 4) - 5.0 * pow(ll, 2)) * sk[3] / 4.0;
    hh[2][1] = -3.0*pow(ll,2)*(-1.0+pow(nn,2)+pow(ll,2))*sqrt(5.0)*nn*sk[1]+(6.0*pow(ll,2)*pow(nn,2)-pow(nn,2)+1.0+6.0*pow(ll,4)-6.0*pow(ll,2))
              *sqrt(10.0)*nn*sk[2]/2.0-(nn*(3.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+1.0-3.0*pow(ll,2)+3.0*pow(ll,4))*sk[3]);
    hh[3][1] = -3.0/4.0*ll*(-1.0+pow(nn,2)+pow(ll,2))*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]+((30.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)
              -27.0*pow(nn,2)-2.0*pow(ll,2)+1.0)*ll*sk[2])/4.0-ll*sqrt(10.0)*(3.0*pow(nn,4)+3.0*pow(ll,2)*pow(nn,2)+pow(ll,2)-1.0)*sk[3]/4.0;
    hh[4][1] = ll*mm*sqrt(3.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1]/2.0-(5.0*pow(nn,2)-1.0)*nn*ll*mm*sqrt(6.0)*sk[2]/2.0
              +ll*mm*(pow(nn,2)+1.0)*sqrt(15.0)*nn*sk[3]/2.0;
    hh[5][1] = 3.0/4.0*(pow(ll,2))*mm*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]-(30.0*pow(ll,2)*pow(nn,2)-5.0*pow(nn,2)-2.0*pow(ll,2)+1.0)*mm
              *sk[2]/4.0+mm*sqrt(10.0)*(3.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+pow(ll,2))*sk[3]/4.0;
    hh[6][1] = 3.0/2.0*ll*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(5.0)*nn*sk[1]-3.0/2.0*nn*sqrt(10.0)*mm*ll
              *(2.0*pow(ll,2)-1.0+pow(nn,2))*sk[2]+3.0/2.0*ll*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*nn*sk[3];
    hh[7][1] = (pow(ll,2))*mm*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*sqrt(5.0)*sk[1]/ 4.0-sqrt(15.0)*mm
              *(6.0*pow(ll,2)*pow(nn,2)-pow(nn,2)+1.0+8.0*pow(ll,4)-6.0*pow(ll,2))
              *sk[2]/4.0+sqrt(6.0)*mm*(3.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+4.0*pow(ll,4)-3.0*pow(ll,2))*sk[3]/4.0;
    hh[1][2] = -(-1.0+pow(nn,2)+pow(ll,2))*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0+sqrt(15.0)*nn*(2.0*pow(nn,4)
              -3.0*pow(nn,2)+10.0*pow(ll,2)*pow(nn,2)+1.0+8.0*pow(ll,4)-8.0*pow(ll,2))*sk[2]/4.0-sqrt(6.0)*nn*(pow(nn,4)+5.0*pow(ll,2)
              *pow(nn,2)-1.0+4.0*pow(ll,4)-pow(ll,2))*sk[3]/4.0;
    hh[2][2] = -3.0*(-1.0+pow(nn,2)+pow(ll,2))*ll*(pow(nn,2))*sqrt(5.0)*sk[1]+(6.0*pow(nn,4)-6.0*pow(nn,2)+6.0*pow(ll,2)*pow(nn,2)+1.0-pow(ll,2))*ll
              *sqrt(10.0)*sk[2]/2.0-(ll*(3.0*pow(nn,4)+3.0*pow(ll,2)*pow(nn,2)-3.0*pow(nn,2)+1.0-2.0*pow(ll,2))*sk[3]);
    hh[3][2] = -3.0/4.0*(-1.0+pow(nn,2)+pow(ll,2))*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]+((30.0*pow(nn,4)-37.0*pow(nn,2)
              +30.0*pow(ll,2)*pow(nn,2)+11.0-12.0*pow(ll,2))*nn*sk[2])/4.0-(-1.0+nn)*sqrt(10.0)*(3.0*pow(nn,3)+3.0*pow(nn,2)-nn
              +3.0*pow(ll,2)*nn-1.0+3.0*pow(ll,2))*nn*sk[3]/4.0;
    hh[4][2] = mm*sqrt(3.0)*(pow(nn,2))*(5.0*pow(nn,2)-3.0)*sk[1]/2.0-(5.0*pow(nn,2)-1.0)*sqrt(3.0)*sqrt(2.0)
              *(2.0*pow(nn,2)-1.0)*mm*sk[2]/4.0+(-1.0+pow(nn,2))*sqrt(15.0)*(pow(nn,2))*mm*sk[3]/2.0;
    hh[5][2] = 3.0/4.0*mm*ll*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]-3.0/2.0*(5.0*pow(nn,2)-2.0)*mm*ll*nn
              *sk[2]+3.0/4.0*nn*sqrt(10.0)*ll*mm*(-1.0+pow(nn,2))*sk[3];
    hh[6][2] = 3.0/2.0*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2))*sqrt(5.0)*sk[1]-(6.0*pow(nn,4)+12.0*pow(ll,2)*pow(nn,2)-5.0*pow(nn,2)
              +1.0-2.0*pow(ll,2))*mm*sqrt(10.0)*sk[2]/4.0+mm*(3.0*pow(nn,4)-pow(nn,2)+6.0*pow(ll,2)*pow(nn,2)-4.0*pow(ll,2))*sk[3]/2.0;
    hh[7][2] = mm*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0-mm*ll*sqrt(15.0)*nn
              *(3.0*pow(nn,2)-2.0+4.0*pow(ll,2))*sk[2]/2.0+sqrt(6.0)*mm*ll*nn*(3.0*pow(nn,2)+1.0+4.0*pow(ll,2))*sk[3]/4.0;
    hh[1][3] = sqrt(2.0)*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*(3.0*pow(nn,2)-1.0)*sqrt(5.0)*sk[1]/8.0-3.0/4.0*(pow(nn,2))*mm
              *(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(5.0)*sk[2]+3.0/8.0*sqrt(2.0)*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2)+1.0)*sk[3];
    hh[2][3] = ll*mm*(3.0*pow(nn,2)-1.0)*sqrt(15.0)*nn*sk[1]/2.0-(3.0*pow(nn,2)-1.0)*ll*nn*mm*sqrt(30.0)*sk[2]/2.0
              +sqrt(3.0)*ll*mm*nn*(3.0*pow(nn,2)-1.0)*sk[3]/2.0;
    hh[3][3] = sqrt(2.0)*mm*(3.0*pow(nn,2)-1.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sk[1]/8.0-(15.0*pow(nn,2)-11.0)*mm*(pow(nn,2))*sqrt(3.0)
              *sk[2]/4.0+(3.0*pow(nn,3)+3.0*pow(nn,2)-nn-1.0)*(-1.0+nn)*sqrt(5.0)*sqrt(2.0)*mm*sqrt(3.0)*sk[3]/8.0;
    hh[4][3] = ((3.0*pow(nn,2)-1.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1])/4.0-3.0/4.0*(5.0*pow(nn,2)-1.0)*(-1.0+pow(nn,2))*nn
              *sqrt(2.0)*sk[2]+3.0/4.0*(pow(-1.0+pow(nn,2),2))*sqrt(5.0)*nn*sk[3];
    hh[5][3] = sqrt(2.0)*ll*(3.0*pow(nn,2)-1.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sk[1]/8.0-(15.0*pow(nn,2)-11.0)*ll*(pow(nn,2))
              *sqrt(3.0)*sk[2]/4.0+(3.0*pow(nn,3)+3.0*pow(nn,2)-nn-1.0)*(-1.0+nn)*sqrt(5.0)*sqrt(2.0)*ll*sqrt(3.0)*sk[3]/8.0;
    hh[6][3] = (2.0*pow(ll,2)-1.0+pow(nn,2))*(3.0*pow(nn,2)-1.0)*sqrt(15.0)*nn*sk[1]/4.0-(3.0*pow(nn,2)-1.0)*(2.0*pow(ll,2)-1.0+pow(nn,2))
              *nn*sqrt(30.0)*sk[2]/4.0+sqrt(3.0)*(2.0*pow(ll,2)-1.0+pow(nn,2))*nn*(3.0*pow(nn,2)-1.0)*sk[3]/4.0;
    hh[7][3] = sqrt(2.0)*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*(3.0*pow(nn,2)-1.0)*sqrt(5.0)*sk[1]/8.0-3.0/4.0*(pow(nn,2))*ll
              *(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(5.0)*sk[2]+ 3.0/8.0*sqrt(2.0)*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*(pow(nn,2)+1.0)*sk[3];
    hh[1][4] = ll*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0-ll*mm*sqrt(15.0)*nn
              *(pow(nn,2)-2.0+4.0*pow(ll,2))*sk[2]/2.0+sqrt(6.0)*mm*ll*nn*(pow(nn,2)+4.0*pow(ll,2)-5.0)*sk[3]/4.0;
    hh[2][4] = 3.0*(pow(ll,2))*mm*(pow(nn,2))*sqrt(5.0)*sk[1]-(6.0*pow(ll,2.0)*pow(nn,2)-pow(nn,2)-pow(ll,2))*mm*sqrt(10.0)*sk[2]
              /2.0+mm*(-2.0*pow(nn,2)+3.0*pow(ll,2.0)*pow(nn,2)+1.0-2.0*pow(ll,2))*sk[3];
    hh[3][4] = 3.0/4.0*ll*mm*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]-3.0/2.0*(5.0*pow(nn,2)-2.0)*ll*mm*nn
              *sk[2]+3.0/4.0*nn*sqrt(10.0)*mm*ll*(-1.0+pow(nn,2))*sk[3];
    hh[4][4] = ll*sqrt(3.0)*(pow(nn,2))*(5.0*pow(nn,2)-3.0)*sk[1]/2.0-(5.0*pow(nn,2)-1.0)*sqrt(3.0)*sqrt(2.0)*(2.0*pow(nn,2)-1.0)
              *ll*sk[2]/4.0+(-1.0+pow(nn,2))*sqrt(15.0)*(pow(nn,2))*ll*sk[3]/2.0;
    hh[5][4] = 3.0/ 4.0*pow(ll,2.0)*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]-(30.0*pow(ll,2.0)*(pow(nn,2))-(5.0*pow(nn,2))
              -12.0*pow(ll,2)+1.0)*nn*sk[2]/4.0+(-1.0+nn)*sqrt(10.0)*(-(2.0*nn)+3.0*pow(ll,2.0)*nn-2.0+3.0*pow(ll,2))*nn*sk[3]/4.0;
    hh[6][4] = 3.0/2.0*ll*(2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2))*sqrt(5.0)*sk[1]-(6.0*pow(nn,4)+12.0*pow(ll,2.0)*pow(nn,2)
              -9.0*pow(nn,2)+1.0-2.0*pow(ll,2))*ll*sqrt(10.0)*sk[2]/4.0+(ll*(3.0*pow(nn,4)-9.0*pow(nn,2)+6.0*pow(ll,2.0)*pow(nn,2)
              +4.0-4.0*pow(ll,2))*sk[3])/2.0;
    hh[7][4] = (pow(ll,2))*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0-sqrt(15.0)*nn*(6.0*pow(ll,2.0)*pow(nn,2)
              -pow(nn,2)+1.0+8.0*pow(ll,4)-8.0*pow(ll,2))*sk[2]/4.0+sqrt(6.0)*nn*(3.0*pow(ll,2.0)*pow(nn,2)-2.0*pow(nn,2)-7.0*pow(ll,2)
              +2.0+4.0*pow(ll,4))*sk[3]/4.0;
    hh[1][5] = (2.0*pow(ll,2)-1.0+pow(nn,2))*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*sqrt(5.0)*sk[1]/8.0-sqrt(15.0)*mm
              *(pow(nn,4)-pow(nn,2)+6.0*pow(ll,2.0)*pow(nn,2)+8.0*pow(ll,4)-6.0*pow(ll,2))*sk[2]/4.0+sqrt(6.0)*mm*(pow(nn,4)
              +6.0*pow(ll,2.0)*pow(nn,2)+2.0*pow(nn,2)+1.0+8.0*pow(ll,4)-6.0*pow(ll,2))*sk[3]/8.0;
    hh[2][5] = 3.0/2.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*mm*sqrt(5.0)*nn*sk[1]-3.0/2.0*nn*sqrt(10.0)*(2.0*pow(ll,2)-1.0
              +pow(nn,2))*ll*mm*sk[2]+ 3.0/2.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*mm*nn*sk[3];
    hh[3][5] = 3.0/8.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*mm*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]-(15.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)
              -11.0*pow(nn,2)-2.0*pow(ll,2))*mm*sk[2]/4.0+mm*sqrt(10.0)*(3.0*pow(nn,4)+2.0*pow(nn,2)+6.0*pow(ll,2)*pow(nn,2)-1.0
              +2.0*pow(ll,2))*sk[3]/8.0;
    hh[4][5] = (2.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(3.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1]/4.0-(5.0*pow(nn,2)-1.0)*nn*(2.0*pow(ll,2)-1.0
              +pow(nn,2))*sqrt(6.0)*sk[2]/4.0+(2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2)+1.0)*sqrt(15.0)*nn*sk[3]/4.0;
    hh[5][5] = 3.0/8.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]-((15.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)
              -21.0*pow(nn,2)+2.0-2.0*pow(ll,2))*ll*sk[2])/4.0+ll*sqrt(10.0)*(3.0*pow(nn,4)+6.0*pow(ll,2)*pow(nn,2)-6.0*pow(nn,2)
              +2.0*pow(ll,2)-1.0)*sk[3]/8.0;
    hh[6][5] = 3.0/4.0*(pow(2.0*pow(ll,2)-1.0+pow(nn,2),2))*sqrt(5.0)*nn*sk[1]-(3.0*pow(nn,4)+12.0*pow(ll,2)*pow(nn,2)-4.0*pow(nn,2)
              +12.0*pow(ll,4)+1.0-12.0*pow(ll,2))*sqrt(10.0)*nn*sk[2]/4.0+(nn*(3.0*pow(nn,4)+12.0*pow(ll,2)*pow(nn,2)+2.0*pow(nn,2)
              -12.0*pow(ll,2)-1.0+12.0*pow(ll,4))*sk[3])/4.0;
    hh[7][5] = (2.0*pow(ll,2)-1.0+pow(nn,2))*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*sqrt(5.0)*sk[1]/8.0-sqrt(15.0)*ll
              *(3.0*pow(nn,4)+10.0*pow(ll,2)*pow(nn,2)-5.0*pow(nn,2)-10.0*pow(ll,2)+8.0*pow(ll,4)+2.0)*sk[2]/4.0+sqrt(6.0)*ll
              *(3.0*pow(nn,4)+10.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+3.0+8.0*pow(ll,4)-10.0*pow(ll,2))*sk[3]/8.0;
  }

  // rotation routine for interaction of an f orbital with an f orbital
  void df(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 3);
    assert(hh.size() >= 5 && hh[0].size() >= 7);

    hh[1][1] = -ll*(-1.0+pow(nn,2)+pow(ll,2))*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*sqrt(5.0)*sk[1]/4.0+sqrt(15.0)*ll
              *(2.0*pow(nn,4)-5.0*pow(nn,2)+10.0*pow(ll,2)*pow(nn,2)+3.0-10.0*pow(ll,2)+8.0*pow(ll,4))*sk[2]/4.0-sqrt(6.0)*ll
              *(pow(nn,4)+5.0*pow(ll,2)*pow(nn,2)-4.0*pow(nn,2)+1.0+4.0*pow(ll,4)-5.0*pow(ll,2))*sk[3]/4.0;
    hh[2][1] = -3.0*(pow(ll,2))*(-1.0+pow(nn,2)+pow(ll,2))*sqrt(5.0)*nn*sk[1]+(6.0*pow(ll,2)*pow(nn,2)-pow(nn,2)+1.0+6.0*pow(ll,4)
              -6.0*pow(ll,2))*sqrt(10.0)*nn*sk[2]/2.0-(nn*(3.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+1.0-3.0*pow(ll,2)+3.0*pow(ll,4))*sk[3]);
    hh[3][1] = -3.0/4.0*ll*(-1.0+pow(nn,2)+pow(ll,2))*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]+((30.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)
              -27.0*pow(nn,2)-2.0*pow(ll,2)+1.0)*ll*sk[2])/4.0-ll*sqrt(10.0)*(3.0*pow(nn,4)+3.0*pow(ll,2)*pow(nn,2)+pow(ll,2)-1.0)*sk[3]/4.0;
    hh[4][1] = ll*mm*sqrt(3.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1]/2.0-(5.0*pow(nn,2)-1.0)*nn*ll*mm*sqrt(6.0)*sk[2]/2.0
              +ll*mm*(pow(nn,2)+1.0)*sqrt(15.0)*nn*sk[3]/2.0;
    hh[5][1] = 3.0/4.0*(pow(ll,2))*mm*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]-(30.0*pow(ll,2)*pow(nn,2)-5.0*pow(nn,2)-2.0*pow(ll,2)+1.0)
              *mm*sk[2]/4.0+mm*sqrt(10.0)*(3.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+pow(ll,2))*sk[3]/4.0;
    hh[6][1] = 3.0/2.0*ll*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(5.0)*nn*sk[1]-3.0/2.0*nn*sqrt(10.0)*mm*ll*(2.0*pow(ll,2)-1.0+pow(nn,2))
              *sk[2]+3.0/2.0*ll*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*nn*sk[3];
    hh[7][1] = (pow(ll,2))*mm*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*sqrt(5.0)*sk[1]/ 4.0-sqrt(15.0)*mm*(6.0*pow(ll,2)*pow(nn,2)
              -pow(nn,2)+1.0+8.0*pow(ll,4)-6.0*pow(ll,2))*sk[2]/4.0+sqrt(6.0)*mm*(3.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+4.0*pow(ll,4)
              -3.0*pow(ll,2))*sk[3]/4.0;
    hh[1][2] = -(-1.0+pow(nn,2)+pow(ll,2))*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0+sqrt(15.0)*nn
              *(2.0*pow(nn,4)-3.0*pow(nn,2)+10.0*pow(ll,2)*pow(nn,2)+1.0+8.0*pow(ll,4)-8.0*pow(ll,2))*sk[2]/4.0-sqrt(6.0)*nn
              *(pow(nn,4)+5.0*pow(ll,2)*pow(nn,2)-1.0+4.0*pow(ll,4)-pow(ll,2))*sk[3]/4.0;
    hh[2][2] = -3.0*(-1.0+pow(nn,2)+pow(ll,2))*ll*(pow(nn,2))*sqrt(5.0)*sk[1]+(6.0*pow(nn,4)-6.0*pow(nn,2)+6.0*pow(ll,2)*pow(nn,2)+1.0
              -pow(ll,2))*ll*sqrt(10.0)*sk[2]/2.0-(ll*(3.0*pow(nn,4)+3.0*pow(ll,2)*pow(nn,2)-3.0*pow(nn,2)+1.0-2.0*pow(ll,2))*sk[3]);
    hh[3][2] = -3.0/4.0*(-1.0+pow(nn,2)+pow(ll,2))*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]+((30.0*pow(nn,4)-37.0*pow(nn,2)
              +30.0*pow(ll,2)*pow(nn,2)+11.0-12.0*pow(ll,2))*nn*sk[2])/4.0-(-1.0+nn)*sqrt(10.0)*(3.0*pow(nn,3)+3.0*pow(nn,2)-nn
              +3.0*pow(ll,2)*nn-1.0+3.0*pow(ll,2))*nn*sk[3]/4.0;
    hh[4][2] = mm*sqrt(3.0)*(pow(nn,2))*(5.0*pow(nn,2)-3.0)*sk[1]/2.0-(5.0*pow(nn,2)-1.0)*sqrt(3.0)*sqrt(2.0)*(2.0*pow(nn,2)-1.0)
              *mm*sk[2]/4.0+(-1.0+pow(nn,2))*sqrt(15.0)*(pow(nn,2))*mm*sk[3]/2.0;
    hh[5][2] = 3.0/4.0*mm*ll*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]-3.0/2.0*(5.0*pow(nn,2)-2.0)*mm*ll*nn*sk[2]
              +3.0/4.0*nn*sqrt(10.0)*ll*mm*(-1.0+pow(nn,2))*sk[3];
    hh[6][2] = 3.0/2.0*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2))*sqrt(5.0)*sk[1]-(6.0*pow(nn,4)+12.0*pow(ll,2)*pow(nn,2)
              -5.0*pow(nn,2)+1.0-2.0*pow(ll,2))*mm*sqrt(10.0)*sk[2]/4.0+mm*(3.0*pow(nn,4)-pow(nn,2)+6.0*pow(ll,2)*pow(nn,2)
              -4.0*pow(ll,2))*sk[3]/2.0;
    hh[7][2] = mm*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0-mm*ll*sqrt(15.0)*nn*(3.0*pow(nn,2)-2.0
              +4.0*pow(ll,2))*sk[2]/2.0+sqrt(6.0)*mm*ll*nn*(3.0*pow(nn,2)+1.0+4.0*pow(ll,2))*sk[3]/4.0;
    hh[1][3] = sqrt(2.0)*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*(3.0*pow(nn,2)-1.0)*sqrt(5.0)*sk[1]/8.0-3.0/4.0*(pow(nn,2))*mm
              *(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(5.0)*sk[2]+3.0/8.0*sqrt(2.0)*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2)+1.0)*sk[3];
    hh[2][3] = ll*mm*(3.0*pow(nn,2)-1.0)*sqrt(15.0)*nn*sk[1]/2.0-(3.0*pow(nn,2)-1.0)*ll*nn*mm*sqrt(30.0)*sk[2]/2.0
              +sqrt(3.0)*ll*mm*nn*(3.0*pow(nn,2)-1.0)*sk[3]/2.0;
    hh[3][3] = sqrt(2.0)*mm*(3.0*pow(nn,2)-1.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sk[1]/8.0-(15.0*pow(nn,2)-11.0)*mm*(pow(nn,2))
              *sqrt(3.0)*sk[2]/4.0+(3.0*pow(nn,3)+3.0*pow(nn,2)-nn-1.0)*(-1.0+nn)*sqrt(5.0)*sqrt(2.0)*mm*sqrt(3.0)*sk[3]/8.0;
    hh[4][3] = ((3.0*pow(nn,2)-1.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1])/4.0-3.0/4.0*(5.0*pow(nn,2)-1.0)*(-1.0+pow(nn,2))*nn*sqrt(2.0)*sk[2]
              +3.0/4.0*pow(-1.0+pow(nn,2),2)*sqrt(5.0)*nn*sk[3];
    hh[5][3] = sqrt(2.0)*ll*(3.0*pow(nn,2)-1.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sk[1]/8.0-(15.0*pow(nn,2)-11.0)*ll*(pow(nn,2))*sqrt(3.0)
              *sk[2]/4.0+(3.0*pow(nn,3)+3.0*pow(nn,2)-nn-1.0)*(-1.0+nn)*sqrt(5.0)*sqrt(2.0)*ll*sqrt(3.0)*sk[3]/8.0;
    hh[6][3] = (2.0*pow(ll,2)-1.0+pow(nn,2))*(3.0*pow(nn,2)-1.0)*sqrt(15.0)*nn*sk[1]/4.0-(3.0*pow(nn,2)-1.0)*(2.0*pow(ll,2)-1.0
              +pow(nn,2))*nn*sqrt(30.0)*sk[2]/4.0+sqrt(3.0)*(2.0*pow(ll,2)-1.0+pow(nn,2))*nn*(3.0*pow(nn,2)-1.0)*sk[3]/4.0;
    hh[7][3] = sqrt(2.0)*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*(3.0*pow(nn,2)-1.0)*sqrt(5.0)*sk[1]/8.0-3.0/4.0*(pow(nn,2))*ll
              *(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(5.0)*sk[2]+ 3.0/8.0*sqrt(2.0)*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*(pow(nn,2)+1.0)*sk[3];
    hh[1][4] = ll*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0-ll*mm*sqrt(15.0)*nn*(pow(nn,2)-2.0+4.0*pow(ll,2))
              *sk[2]/2.0+sqrt(6.0)*mm*ll*nn*(pow(nn,2)+4.0*pow(ll,2)-5.0)*sk[3]/4.0;
    hh[2][4] = 3.0*(pow(ll,2))*mm*(pow(nn,2))*sqrt(5.0)*sk[1]-(6.0*pow(ll,2.0)*pow(nn,2)-pow(nn,2)-pow(ll,2))*mm*sqrt(10.0)*sk[2]
              /2.0+mm*(-2.0*pow(nn,2)+3.0*pow(ll,2.0)*pow(nn,2)+1.0-2.0*pow(ll,2))*sk[3];
    hh[3][4] = 3.0/4.0*ll*mm*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]-3.0/2.0*(5.0*pow(nn,2)-2.0)*ll*mm*nn*sk[2]+3.0/4.0*nn*sqrt(10.0)
              *mm*ll*(-1.0+pow(nn,2))*sk[3];
    hh[4][4] = ll*sqrt(3.0)*(pow(nn,2))*(5.0*pow(nn,2)-3.0)*sk[1]/2.0-(5.0*pow(nn,2)-1.0)*sqrt(3.0)*sqrt(2.0)*(2.0*pow(nn,2)-1.0)
              *ll*sk[2]/4.0+(-1.0+pow(nn,2))*sqrt(15.0)*(pow(nn,2))*ll*sk[3]/2.0;
    hh[5][4] = 3.0/ 4.0*pow(ll,2.0)*sqrt(2.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]-(30.0*pow(ll,2.0)*(pow(nn,2))-(5.0*pow(nn,2))
              -12.0*pow(ll,2)+1.0)*nn*sk[2]/4.0+(-1.0+nn)*sqrt(10.0)*(-(2.0*nn)+3.0*pow(ll,2.0)*nn-2.0+3.0*pow(ll,2))*nn*sk[3]/4.0;
    hh[6][4] = 3.0/2.0*ll*(2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2))*sqrt(5.0)*sk[1]-(6.0*pow(nn,4)+12.0*pow(ll,2.0)*pow(nn,2)
              -9.0*pow(nn,2)+1.0-2.0*pow(ll,2))*ll*sqrt(10.0)*sk[2]/4.0+(ll*(3.0*pow(nn,4)-9.0*pow(nn,2)+6.0*pow(ll,2.0)*pow(nn,2)
              +4.0-4.0*pow(ll,2))*sk[3])/2.0;
    hh[7][4] = (pow(ll,2))*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*nn*sqrt(5.0)*sk[1]/4.0-sqrt(15.0)*nn*(6.0*pow(ll,2.0)*pow(nn,2)
              -pow(nn,2)+1.0+8.0*pow(ll,4)-8.0*pow(ll,2))*sk[2]/4.0+sqrt(6.0)*nn*(3.0*pow(ll,2.0)*pow(nn,2)-2.0*pow(nn,2)-7.0*pow(ll,2)
              +2.0+4.0*pow(ll,4))*sk[3]/4.0;
    hh[1][5] = (2.0*pow(ll,2)-1.0+pow(nn,2))*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*sqrt(5.0)*sk[1]/8.0-sqrt(15.0)*mm
              *(pow(nn,4)-pow(nn,2)+6.0*pow(ll,2.0)*pow(nn,2)+8.0*pow(ll,4)-6.0*pow(ll,2))*sk[2]/4.0+sqrt(6.0)*mm*(pow(nn,4)
              +6.0*pow(ll,2.0)*pow(nn,2)+2.0*pow(nn,2)+1.0+8.0*pow(ll,4)-6.0*pow(ll,2))*sk[3]/8.0;
    hh[2][5] = 3.0/2.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*mm*sqrt(5.0)*nn*sk[1]-3.0/2.0*nn*sqrt(10.0)*(2.0*pow(ll,2)-1.0
              +pow(nn,2))*ll*mm*sk[2]+ 3.0/2.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*mm*nn*sk[3];
    hh[3][5] = 3.0/8.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*mm*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]-(15.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)
              -11.0*pow(nn,2)-2.0*pow(ll,2))*mm*sk[2]/4.0+mm*sqrt(10.0)*(3.0*pow(nn,4)+2.0*pow(nn,2)+6.0*pow(ll,2)*pow(nn,2)-1.0
              +2.0*pow(ll,2))*sk[3]/8.0;
    hh[4][5] = (2.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(3.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1]/4.0-(5.0*pow(nn,2)-1.0)*nn*(2.0*pow(ll,2)-1.0
              +pow(nn,2))*sqrt(6.0)*sk[2]/4.0+(2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2)+1.0)*sqrt(15.0)*nn*sk[3]/4.0;
    hh[5][5] = 3.0/8.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[1]-((15.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)
              -21.0*pow(nn,2)+2.0-2.0*pow(ll,2))*ll*sk[2])/4.0+ll*sqrt(10.0)*(3.0*pow(nn,4)+6.0*pow(ll,2)*pow(nn,2)-6.0*pow(nn,2)
              +2.0*pow(ll,2)-1.0)*sk[3]/8.0;
    hh[6][5] = 3.0/4.0*pow((2.0*pow(ll,2)-1.0+pow(nn,2)),2)*sqrt(5.0)*nn*sk[1]-(3.0*pow(nn,4)+12.0*pow(ll,2)*pow(nn,2)-4.0*pow(nn,2)
              +12.0*pow(ll,4)+1.0-12.0*pow(ll,2))*sqrt(10.0)*nn*sk[2]/4.0+(nn*(3.0*pow(nn,4)+12.0*pow(ll,2)*pow(nn,2)+2.0*pow(nn,2)
              -12.0*pow(ll,2)-1.0+12.0*pow(ll,4))*sk[3])/4.0;
    hh[7][5] = (2.0*pow(ll,2)-1.0+pow(nn,2))*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*sqrt(5.0)*sk[1]/8.0-sqrt(15.0)*ll
              *(3.0*pow(nn,4)+10.0*pow(ll,2)*pow(nn,2)-5.0*pow(nn,2)-10.0*pow(ll,2)+8.0*pow(ll,4)+2.0)*sk[2]/4.0+sqrt(6.0)*ll
              *(3.0*pow(nn,4)+10.0*pow(ll,2)*pow(nn,2)-2.0*pow(nn,2)+3.0+8.0*pow(ll,4)-10.0*pow(ll,2))*sk[3]/8.0;
  }

  // rotation routine for interaction of an f orbital with an f orbital
  void ff(

    // dimeric block to put the results in to
    vector<vector<double>> &hh,

    // directional cosine ll
    double const &ll,

    // directional cosine mm
    double const &mm,

    // directional cosine nn
    double const &nn,

    // Slater-Koster table for dimer element of the Slater-Koster table
    vector<double> const &sk

  )
  {
    assert(sk.size() == 4);
    assert(hh.size() >= 7 && hh[0].size() >= 7);

    hh[1][1] = -5.0/8.0*(-1.0+pow(nn,2)+pow(ll,2))*pow(4.0*pow(ll,2)-1.0+pow(nn,2),2)*sk[1]+(15.0/16.0*(pow(nn,6))- 15.0
              /8.0*(pow(nn,4))+135.0/ 16.0*(pow(nn,4))*(pow(ll,2))-135.0/8.0*(pow(ll,2))*(pow(nn,2))+15.0/16.0*(pow(nn,2))
              +45.0/2.0*(pow(nn,2))*(pow(ll,4))+135.0/16.0*(pow(ll,2))-45.0/2.0*(pow(ll,4))+(15.0*pow(ll,6)))*sk[2]+(-3.0/ 8.0
              *(pow(nn,6))-27.0/8.0*(pow(nn,4))*(pow(ll,2))-3.0/8.0*(pow(nn,4))+3.0/8.0*(pow(nn,2))-(9.0*pow(nn,2)*pow(ll,4))
              +27.0/4.0*(pow(ll,2))*(pow(nn,2))+3.0/8.0-(6.0*pow(ll,6))-27.0/8.0*(pow(ll,2))+(9.0*pow(ll,4)))*sk[3]+((pow(nn,6))/16.0
              +3.0/8.0*(pow(nn,4))+9.0/16.0*(pow(nn,4))*(pow(ll,2))-9.0/8.0*(pow(ll,2))*(pow(nn,2))+9.0/16.0*(pow(nn,2))
              +3.0/2.0*(pow(nn,2))*(pow(ll,4))+9.0/16.0*(pow(ll,2))-3.0/2.0*(pow(ll,4))+(pow(ll,6)))*sk[4];
    hh[2][1] = - 5.0/4.0*(-1.0+pow(nn,2)+pow(ll,2))*(4.0*pow(ll,2)-1.0+pow(nn,2))*ll*sqrt(6.0)*nn*sk[1]+5.0/8.0*sqrt(6.0)*nn*ll
              *(3.0*pow(nn,4)+15.0*pow(ll,2)*pow(nn,2)-7.0*pow(nn,2)+4.0-15.0*pow(ll,2)+12.0*pow(ll,4))*sk[2]-sqrt(6.0)*nn*ll
              *(3.0*pow(nn,4)+15.0*pow(ll,2)*pow(nn,2)-10.0*pow(nn,2)+5.0-15.0*pow(ll,2)+12.0*pow(ll,4))*sk[3]/4.0+ll*sqrt(6.0)
              *(pow(nn,4)+5.0*pow(ll,2)*pow(nn,2)-5.0*pow(nn,2)+4.0*pow(ll,4)-5.0*pow(ll,2))*nn*sk[4]/8.0;
    hh[3][1] = -(-1.0+pow(nn,2)+pow(ll,2))*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(5.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sk[1]/8.0
              +sqrt(15.0)*(15.0*pow(nn,6)-26.0*pow(nn,4)+75.0*pow(ll,2)*pow(nn,4)-70.0*pow(ll,2)*pow(nn,2)+11.0*pow(nn,2)
              +60.0*pow(ll,4)*pow(nn,2)-4.0*pow(ll,4)+3.0*pow(ll,2))*sk[2]/16.0-sqrt(15.0)*(3.0*pow(nn,6)-pow(nn,4)
              +15.0*pow(ll,2)*pow(nn,4)-3.0*pow(nn,2)-2.0*pow(ll,2)*pow(nn,2)+12.0*pow(ll,4)*pow(nn,2)+1.0+4.0*pow(ll,4)
              -5.0*pow(ll,2))*sk[3]/8.0+sqrt(15.0)*(pow(nn,6)+2.0*pow(nn,4)+5.0*pow(ll,2)*pow(nn,4)-3.0*pow(nn,2)
              +6.0*pow(ll,2)*pow(nn,2)+4.0*pow(ll,4)*pow(nn,2)-3.0*pow(ll,2)+4.0*pow(ll,4))*sk[4]/16.0;
    hh[4][1] = mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(2.0)*sqrt(5.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1]/8.0-3.0/16.0*mm*(4.0*pow(ll,2)
              -1.0+pow(nn,2))*sqrt(5.0)*nn*sqrt(2.0)*(5.0*pow(nn,2)-1.0)*sk[2]+3.0/8.0*sqrt(10.0)*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))
              *(pow(nn,2)+1.0)*nn*sk[3]-sqrt(5.0)*sqrt(2.0)*(pow(nn,2)+3.0)*nn*(4.0*pow(ll,2)-1.0+pow(nn,2))*mm*sk[4]/16.0;
    hh[5][1] = mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*ll*sqrt(5.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sk[1]/8.0-mm*sqrt(15.0)*ll*(15.0*pow(nn,4)
              +60.0*pow(ll,2)*pow(nn,2)-26.0*pow(nn,2)+3.0-4.0*pow(ll,2))*sk[2]/ 16.0+mm*sqrt(15.0)*ll*(3.0*pow(nn,4)
              +12.0*pow(ll,2)*pow(nn,2)-10.0*pow(nn,2)+4.0*pow(ll,2)-1.0)*sk[3]/8.0-mm*sqrt(15.0)*ll*(pow(nn,4)+4.0*pow(ll,2)*pow(nn,2)
              -6.0*pow(nn,2)+4.0*pow(ll,2)-3.0)*sk[4]/16.0;
    hh[6][1] = 5.0/8.0*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*(2.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(6.0)*nn*sk[1]-5.0/16.0*sqrt(6.0)*mm*nn
              *(3.0*pow(nn,4)+18.0*pow(ll,2)*pow(nn,2)-4.0*pow(nn,2)+24.0*pow(ll,4)+1.0-18.0*pow(ll,2))*sk[2]+sqrt(6.0)*mm*nn
              *(3.0*pow(nn,4)+18.0*pow(ll,2)*pow(nn,2)+2.0*pow(nn,2)+24.0*pow(ll,4)-18.0*pow(ll,2)-1.0)*sk[3]/8.0-mm*sqrt(6.0)
              *(pow(nn,4)+6.0*pow(ll,2)*pow(nn,2)+4.0*pow(nn,2)+3.0-6.0*pow(ll,2)+8.0*pow(ll,4))*nn*sk[4]/16.0;
    hh[7][1] = 5.0/8.0*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sk[1]-15.0/16.0*ll*(4.0*pow(ll,2)
              -3.0+3.0*pow(nn,2))*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sk[2]+3.0/8.0*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*mm
              *(4.0*pow(ll,2)-1.0+pow(nn,2))*sk[3]-ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*mm*(4.0*pow(ll,2)-1.0+pow(nn,2))*sk[4]/16.0;
    hh[1][2] = hh[2][1];
    hh[2][2] = -(15.0*pow(ll,2)*(-1.0+pow(nn,2)+pow(ll,2))*pow(nn,2)*sk[1])+(45.0/2.0*(pow(ll,2))*(pow(nn,4))-5.0/2.0*(pow(nn,4))
              -(25.0*pow(ll,2)*pow(nn,2))+45.0/2.0*(pow(ll,4))*(pow(nn,2))+5.0/2.0*(pow(nn,2))+5.0/2.0*(pow(ll,2))-5.0/2.0*(pow(ll,4)))
              *sk[2]+((-9.0*pow(ll,2)*pow(nn,4)+4.0*pow(nn,4)+13.0*pow(ll,2)*pow(nn,2)-9.0*pow(ll,4)*pow(nn,2)-4.0*pow(nn,2)
              -4.0*pow(ll,2)+4.0*pow(ll,4)+1.0)*sk[3])+(3.0/2.0*(pow(ll,2))*(pow(nn,4))-3.0/2.0*(pow(nn,4))+3.0/2.0*(pow(nn,2))
              +3.0/2.0*(pow(ll,4))*(pow(nn,2))-(3.0*pow(ll,2)*pow(nn,2))+3.0/2.0*(pow(ll,2))-3.0/2.0*(pow(ll,4)))*sk[4];
    hh[3][2] = - 3.0/4.0*ll*(-1.0+pow(nn,2)+pow(ll,2))*sqrt(10.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]+(45.0*pow(nn,4)-53.0*pow(nn,2)
              +45.0*pow(ll,2)*pow(nn,2)+12.0-13.0*pow(ll,2))*nn*ll*sqrt(10.0)*sk[2]/8.0-sqrt(10.0)*ll*(9.0*pow(nn,4)+9.0*pow(ll,2)
              *pow(nn,2)-10.0*pow(nn,2)+3.0-5.0*pow(ll,2))*nn*sk[3]/4.0+3.0/8.0*ll*sqrt(2.0)*(-1.0+nn)*sqrt(5.0)*(pow(nn,3)+pow(nn,2)
              +pow(ll,2)*nn+pow(ll,2))*nn*sk[4];
    hh[4][2] = ll*mm*sqrt(15.0)*(pow(nn,2))*(5.0*pow(nn,2)-3.0)*sk[1]/2.0-(5.0*pow(nn,2)-1.0)*(3.0*pow(nn,2)-1.0)*ll*mm
              *sqrt(15.0)*sk[2]/4.0+ll*mm*(pow(nn,2))*(3.0*pow(nn,2)-1.0)*sqrt(15.0)*sk[3]/2.0-ll*mm*sqrt(3.0)*(pow(nn,3)+pow(nn,2)
              +nn+1.0)*sqrt(5.0)*(-1.0+nn)*sk[4]/4.0;
    hh[5][2] = 3.0/4.0*(pow(ll,2))*mm*sqrt(10.0)*nn*(5.0*pow(nn,2)-1.0)*sk[1]-(45.0*pow(ll,2)*pow(nn,2)-5.0*pow(nn,2)
              -13.0*pow(ll,2)+1.0)*nn*mm*sqrt(10.0)*sk[2]/8.0+sqrt(10.0)*mm*(-4.0*pow(nn,2)+9.0*pow(ll,2)*pow(nn,2)+2.0-5.0*pow(ll,2))
              *nn*sk[3]/4.0-3.0/8.0*mm*sqrt(2.0)*(-1.0+nn)*sqrt(5.0)*(-nn+pow(ll,2)*nn-1.0+pow(ll,2))*nn*sk[4];
    hh[6][2] = 15.0/2.0*ll*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2))*sk[1]-5.0/4.0*(9.0*pow(nn,2)-1.0)*(2.0*pow(ll,2)
              -1.0+pow(nn,2))*ll*mm*sk[2]+(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*mm*(9.0*pow(nn,2)-4.0)*sk[3]/2.0-3.0/4.0
              *(-1.0+pow(nn,2))*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*mm*sk[4];
    hh[7][2] = 5.0/4.0*(pow(ll,2))*mm*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*nn*sk[1]-5.0/8.0*sqrt(6.0)*mm*nn
              *(9.0*pow(nn,2)*pow(ll,2)-pow(nn,2)+1.0+12.0*pow(ll,4)-9.0*pow(ll,2))*sk[2]+sqrt(6.0)*mm*nn*(9.0*pow(nn,2)*pow(ll,2)
              -4.0*pow(nn,2)+12.0*pow(ll,4)-9.0*pow(ll,2)+2.0)*sk[3]/4.0-mm*sqrt(6.0)*(3.0*pow(nn,2)*pow(ll,2)-3.0*pow(nn,2)
              -1.0-3.0*pow(ll,2)+4.0*pow(ll,4))*nn*sk[4]/8.0;
    hh[1][3] = hh[3][1];
    hh[2][3] = hh[3][2];
    hh[3][3] = -3.0/8.0*(-1.0+pow(nn,2)+pow(ll,2))*pow((5.0*pow(nn,2)-1.0),2)*sk[1]+(225.0/16.0*(pow(nn,6))-165.0/8.0*(pow(nn,4))
              +225.0/16.0*(pow(ll,2))*(pow(nn,4))+121.0/16.0*(pow(nn,2))-65.0/8.0*(pow(nn,2))*(pow(ll,2))+(pow(ll,2))/16.0)*sk[2]
              -5.0/8.0*(-1.0+nn)*(9.0*pow(nn,5)+9.0*pow(nn,4)-6.0*pow(nn,3)+9.0*pow(ll,2)*pow(nn,3)-6.0*pow(nn,2)+9.0*pow(nn,2)*pow(ll,2)
              +nn-pow(ll,2)*nn+1.0-pow(ll,2))*sk[3]+15.0/16.0*(-pow(nn,2)+pow(nn,4)+pow(nn,2)*pow(ll,2)-pow(ll,2))*(-1.0+pow(nn,2))*sk[4];
    hh[4][3] = mm*sqrt(2.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*nn*(5.0*pow(nn,2)-3.0)*sk[1]/8.0-(5.0*pow(nn,2)-1.0)*sqrt(3.0)*sqrt(2.0)
              *(15.0*pow(nn,2)-11.0)*nn*mm*sk[2]/ 16.0+5.0/8.0*(-1.0+pow(nn,2))*nn*sqrt(3.0)*(3.0*pow(nn,2)-1.0)*sqrt(2.0)*mm*sk[3]
              -5.0/16.0*(-2.0*pow(nn,2)+1.0+pow(nn,4))*sqrt(2.0)*nn*sqrt(3.0)*mm*sk[4];
    hh[5][3] = 3.0/8.0*mm*ll*pow((5.0*pow(nn,2)-1.0),2)*sk[1]-(225.0*pow(nn,4)-130.0*pow(nn,2)+1.0)*mm*ll*sk[2]/16.0
              +5.0/8.0*(9.0*pow(nn,3)+9.0*pow(nn,2)-nn-1.0)*ll*(-1.0+nn)*mm*sk[3]- 15.0/16.0*pow((-1.0+pow(nn,2)),2)*mm*ll*sk[4];
    hh[6][3] = 3.0/8.0*mm*(2.0*pow(ll,2)-1.0+pow(nn,2))*(5.0*pow(nn,2)-1.0)*sqrt(10.0)*nn*sk[1]-(45.0*pow(nn,4)+90.0*pow(nn,2)
              *pow(ll,2)-48.0*pow(nn,2)+11.0-26.0*pow(ll,2))*nn*mm*sqrt(10.0)*sk[2]/16.0+sqrt(10.0)*mm*(9.0*pow(nn,4)-6.0*pow(nn,2)
              +18.0*pow(nn,2)*pow(ll,2)-10.0*pow(ll,2)+1.0)*nn*sk[3]/8.0-3.0/16.0*mm*(-1.0+nn)*sqrt(5.0)*sqrt(2.0)*(pow(nn,3)
              +pow(nn,2)+2.0*nn*pow(ll,2)+nn+2.0*pow(ll,2)+1.0)*nn*sk[4];
    hh[7][3] = mm*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sqrt(5.0)*sk[1]/8.0-mm*ll*sqrt(15.0)
              *(45.0*pow(nn,4)+60.0*pow(nn,2)*pow(ll,2)-38.0*pow(nn,2)+1.0-4.0*pow(ll,2))*sk[2]/16.0+mm*ll*sqrt(15.0)*(9.0*pow(nn,4)
              +2.0*pow(nn,2)+12.0*pow(nn,2)*pow(ll,2)+4.0*pow(ll,2)-3.0)*sk[3]/8.0-mm*sqrt(15.0)*ll*(3.0*pow(nn,4)+4.0*pow(nn,2)
              *pow(ll,2)+6.0*pow(nn,2)-1.0+4.0*pow(ll,2))*sk[4]/ 16.0;
    hh[1][4] = hh[4][1];
    hh[2][4] = hh[4][2];
    hh[3][4] = hh[4][3];
    hh[4][4] = (pow(nn,2)*pow(5.0*pow(nn,2)-3.0,2)*sk[1])/4.0-3.0/8.0*(-1.0+pow(nn,2))*pow((5.0*pow(nn,2)-1.0),2)*sk[2]
              +15.0/4.0*(pow(nn,2))*pow((-1.0+pow(nn,2)),2)*sk[3]-5.0/8.0*(-1.0+pow(nn,2))*(-2.0*pow(nn,2)+1.0+pow(nn,4))*sk[4];
    hh[5][4] = sqrt(2.0)*ll*nn*(5.0*pow(nn,2)-3.0)*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sk[1]/8.0-(15.0*pow(nn,2)-11.0)
              *nn*ll*(5.0*pow(nn,2)-1.0)*sqrt(3.0)*sqrt(2.0)*sk[2]/16.0+5.0/8.0*(3.0*pow(nn,3)+3.0*pow(nn,2)-nn-1.0)*(-1.0+nn)
              *sqrt(2.0)*ll*nn*sqrt(3.0)*sk[3]-5.0/16.0*nn*sqrt(3.0)*ll*pow((-1.0+pow(nn,2)),2)*sqrt(2.0)*sk[4];
    hh[6][4] = (2.0*pow(ll,2)-1.0+pow(nn,2))*(pow(nn,2))*(5.0*pow(nn,2)-3.0)*sqrt(15.0)*sk[1]/4.0-(3.0*pow(nn,2)-1.0)
              *(2.0*pow(ll,2)-1.0+pow(nn,2))*(5.0*pow(nn,2)-1.0)*sqrt(15.0)*sk[2]/8.0+sqrt(15.0)*(pow(nn,2))*(2.0*pow(ll,2)-1.0
              +pow(nn,2))*(3.0*pow(nn,2)-1.0)*sk[3]/4.0-sqrt(5.0)*(2.0*pow(ll,2)-1.0+pow(nn,2))*sqrt(3.0)*(-1.0+pow(nn,4))*sk[4]/8.0;
    hh[7][4] = sqrt(2.0)*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*nn*(5.0*pow(nn,2)-3.0)*sqrt(5.0)*sk[1]/8.0-3.0/16.0*sqrt(2.0)
              *(5.0*pow(nn,2)-1.0)*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(5.0)*nn*sk[2]+3.0/8.0*sqrt(10.0)*nn*ll*(4.0*pow(ll,2)
              -3.0+3.0*pow(nn,2))*(pow(nn,2)+1.0)*sk[3]-sqrt(2.0)*sqrt(5.0)*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*nn*(pow(nn,2)+3.0)*sk[4]/16.0;
    hh[1][5] = hh[5][1];
    hh[2][5] = hh[5][2];
    hh[3][5] = hh[5][3];
    hh[4][5] = hh[5][4];
    hh[5][5] = 3.0/8.0*(pow(ll,2))*pow((5.0*pow(nn,2)-1.0),2)*sk[1]+(-225.0/ 16.0*(pow(ll,2))*(pow(nn,4))+25.0/16.0*(pow(nn,4))
              +65.0/8.0*(pow(ll,2))*(pow(nn,2))-5.0/8.0*(pow(nn,2))+1.0/16.0-(pow(ll,2))/16.0)*sk[2]+5.0/8.0*(-1.0+nn)
              *(9.0*pow(ll,2)*pow(nn,3)-4.0*pow(nn,3)+9.0*pow(ll,2)*pow(nn,2)-4.0*pow(nn,2)-pow(ll,2)*nn
              -pow(ll,2))*sk[3]- 15.0/16.0*(pow(ll,2)*pow(nn,2)+1.0-pow(nn,2)-pow(ll,2))*(-1.0+pow(nn,2))*sk[4];
    hh[6][5] = 3.0/ 8.0*ll*(2.0*pow(ll,2)-1.0+pow(nn,2))*(5.0*pow(nn,2)-1.0)*sqrt(10.0)*nn*sk[1]-(45.0*pow(nn,4)+90.0*pow(ll,2)
              *pow(nn,2)-68.0*pow(nn,2)+15.0-26.0*pow(ll,2))*nn*ll*sqrt(10.0)*sk[2]/ 16.0+sqrt(10.0)*ll*(9.0*pow(nn,4)-22.0*pow(nn,2)
              +18.0*pow(ll,2)*pow(nn,2)+9.0-10.0*pow(ll,2))*nn*sk[3]/ 8.0-3.0/ 16.0*ll*(-1.0+nn)*sqrt(5.0)*sqrt(2.0)*(pow(nn,3)
              +pow(nn,2)+2.0*pow(ll,2)*nn-3.0*nn+2.0*pow(ll,2)-3.0)*nn*sk[4];
    hh[7][5] = (pow(ll,2))*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(3.0)*(5.0*pow(nn,2)-1.0)*sqrt(5.0)*sk[1]/8.0-sqrt(15.0)
              *(-5.0*pow(nn,4)+45.0*pow(ll,2)*pow(nn,4)-58.0*pow(ll,2)*pow(nn,2)+6.0*pow(nn,2)+60.0*pow(ll,4)*pow(nn,2)-1.0
              +5.0*pow(ll,2)-4.0*pow(ll,4))*sk[2]/16.0+sqrt(15.0)*(-4.0*pow(nn,4)+9.0*pow(ll,2)*pow(nn,4)-14.0*pow(ll,2)*pow(nn,2)
              +12.0*pow(ll,4)*pow(nn,2)+4.0*pow(nn,2)-3.0*pow(ll,2)+4.0*pow(ll,4))*sk[3]/8.0-sqrt(15.0)*(3.0*pow(ll,2.0)*pow(nn,4)
              -3.0*pow(nn,4)+2*pow(nn,2)+4.0*pow(ll,4)*pow(nn,2)-6.0*pow(ll,2)*pow(nn,2)+1.0-5.0*pow(ll,2)+4.0*pow(ll,4))*sk[4]/16.0;
    hh[1][6] = hh[6][1];
    hh[2][6] = hh[6][2];
    hh[3][6] = hh[6][3];
    hh[4][6] = hh[6][4];
    hh[5][6] = hh[6][5];
    hh[6][6] = 15.0/4.0*pow((2.0*pow(ll,2)-1.0+pow(nn,2)),2)*(pow(nn,2))*sk[1]+(- 45.0/8.0*(pow(nn,6))-45.0/2.0*(pow(ll,2))*(pow(nn,4))
              +75.0/8.0*(pow(nn,4))+(25.0*pow(ll,2)*pow(nn,2))-35.0/8.0*(pow(nn,2))-45.0/2.0*(pow(ll,4))*(pow(nn,2))-5.0/2.0*(pow(ll,2))
              +5.0/8.0+5.0/2.0*(pow(ll,4)))*sk[2]+(9.0/4.0*(pow(nn,6))+(9.0*pow(ll,2)*pow(nn,4))-3.0/2.0*(pow(nn,4))-(13.0*pow(ll,2)
              *pow(nn,2))+(pow(nn,2))/4.0+(9.0*pow(ll,4)*pow(nn,2))-(4.0*pow(ll,4))+(4.0*pow(ll,2)))*sk[3]+(- 3.0/8.0*(pow(nn,6))
              -3.0/2.0*(pow(ll,2))*(pow(nn,4))-3.0/8.0*(pow(nn,4))+ 3.0/8.0*(pow(nn,2))+(3.0*pow(ll,2)*pow(nn,2))-3.0/2.0*(pow(ll,4))
              *(pow(nn,2))+3.0/2.0*(pow(ll,4))-3.0/2.0*(pow(ll,2))+3.0/8.0)*sk[4];
    hh[7][6] = 5.0/8.0*(2.0*pow(ll,2)-1.0+pow(nn,2))*ll*(4.0*pow(ll,2)-3.0+3.0*pow(nn,2))*sqrt(6.0)*nn*sk[1]-5.0/ 16.0
              *sqrt(6.0)*ll*nn*(9.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)-16.0*pow(nn,2)-30.0*pow(ll,2)+24.0*pow(ll,4)+7.0)*sk[2]
              +sqrt(6.0)*ll*nn*(9.0*pow(nn,4)+30.0*pow(ll,2)*pow(nn,2)-10.0*pow(nn,2)+24.0*pow(ll,4)-30.0*pow(ll,2)+5.0)
              *sk[3]/8.0-sqrt(6.0)*ll*(3.0*pow(nn,4)+10.0*pow(ll,2)*pow(nn,2)+8.0*pow(ll,4)+5.0-10.0*pow(ll,2))*nn*sk[4]/16.0;
    hh[1][7] = hh[7][1];
    hh[2][7] = hh[7][2];
    hh[3][7] = hh[7][3];
    hh[4][7] = hh[7][4];
    hh[5][7] = hh[7][5];
    hh[6][7] = hh[7][6];
    hh[7][7] = 5.0/8.0*(pow(ll,2))*pow((4.0*pow(ll,2)-3.0+3.0*pow(nn,2)),2)*sk[1]+(-135.0/16.0*(pow(ll,2))*(pow(nn,4))+15.0/16.0
              *(pow(nn,4))-45.0/2.0*(pow(ll,4))*(pow(nn,2))-15.0/8.0*(pow(nn,2))+135.0/8.0*(pow(ll,2))*(pow(nn,2))+45.0/2.0*(pow(ll,4))
              +15.0/16.0-135.0/16.0*(pow(ll,2))-(15.0*pow(ll,6)))*sk[2]+(27.0/8.0*(pow(ll,2))*(pow(nn,4))-3.0/2.0*(pow(nn,4))+3.0/2.0
              *(pow(nn,2))+(9.0*pow(ll,4)*pow(nn,2))-27.0/4.0*(pow(ll,2))*(pow(nn,2))+27.0/8.0*(pow(ll,2))-(9.0*pow(ll,4))
              +(6.0*pow(ll,6)))*sk[3]+(9.0/16.0*(pow(nn,4))-9.0/ 16.0*(pow(ll,2))*(pow(nn,4))-3.0/2.0*(pow(ll,4))*(pow(nn,2))
              +3.0/8.0*(pow(nn,2))+9.0/8.0*(pow(ll,2))*(pow(nn,2))+3.0/2.0*(pow(ll,4))+1.0/16.0-9.0/16.0*(pow(ll,2))-(pow(ll,6)))*sk[4];
  }

};

#endif
