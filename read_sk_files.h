#ifndef READ_SK_FILES_H
#define READ_SK_FILES_H

#include <iostream>
#include "t_orbitals.h"
#include "t_slater.h"
#include <vector>
#include <assert.h>
using namespace std;

// Reads Slater-Koster files
// Should be replaced with a more sophisticated routine, once the new SK-format 
// has been established

class TListCharLc{};
//class TSlater{};
class TListIntR1{};
class TRangeSepSKTag{};
double const maxL = 0.0;
// class TOldSKData{};
class TSlakoEqGrid{}; //
class TSplineRepInp{};
class TPolyRepInp{};

// Represents the Slater-Koster data in an SK file.
struct TOldSKData{

  // Grid separation
  double dist;

  // Nr. of grid points
  int nGrid;

  //Atomic eigenvalues.
  double skSelf(); //4???

  // Hubbard Us
  double skHubbU();

  // Occupations
  double skOcc();

  //Mass of the atom
  double mass;

  // Table for H
  vector<vector<double>> skHam;

  // Table for S
  vector<vector<double>> skOver;

};


TOrbitals orb;

void readSKFile(
  
  // List of SK file names to read in for every interaction
  vector<vector<TListCharLc>> &skFiles,

  // Nr. of species in the system
  int const &nSpecies,

  // Data type for slako information
  TSlater &slako,

  // Information about the orbitals in the system
  TOrbitals const &in,

  // For every species, a list of rank one arrays. Each array contains the angular momenta to pick
  // from the appropriate SK-files.
  vector<TListIntR1> &angShells,
  
  // Are the Hubbard Us different for each l-shell?
  bool const &orbRes,

  // Method of the sk interpolation
  int const &skInterMeth,

  // is this a polynomial or spline repulsive?
  vector<vector<bool>> &repPoly,

  // Distances to artificially truncated tables of SK integrals
  double const &truncationCutOff,

  // if calculation range separated then read omega from end of SK file
  TRangeSepSKTag &rangeSepSK

)
{
  int iSp1, iSp2, nSK1, nSK2, iSK1, iSK2, ind, nInteract, iSh1;
  int angShell(maxL+1), nShell;
  bool readRep, readAtomic;
  string filename; //character(lc) :: fileName
  vector<vector<double>> skHam, skOver; // real(dp), allocatable, target :: skHam(:,:), skOver(:,:)
  double dist;
  vector<vector<TOldSKData>> skData12, skData21;
  TSlakoEqGrid pSklakoEqGrid1, pSlakoEqGrid2;
  TSplineRepInp repSplineIn1, repSplineIn2;
  TPolyRepInp repPolyIn1, repPolyIn2;

  // if artificially cutting the SK tables
  int nEntries;

  assert(skFiles.size() == skFiles[0].size());
  assert(skFiles.size() > 0 && skFiles.size() == nSpecies);
  assert(repPoly.size() == skFiles.size() && repPoly[0].size() == skFiles[0].size());

  //slako.skSelf(orb.mShell, nSpecies); No need to allocate space for vector?
  //slako.skHubbU(orb.mShell, nSpecies);
  //slako.skOcc(orb.mShell, nSpecies);
  //slako.mass(nSpecies);
  fill(slako.skSelf.begin(), slako.skSelf.end(), 0.0);
  fill(slako.skHubbU.begin(), slako.skHubbU.end(), 0.0);
  fill(slako.skOcc.begin(), slako.skOcc.end(), 0.0);

  //slako.skHamCont;
  slako.skHamCont(nSpecies);
  //slako.skOverCont;
  slako.skOverCont(nSpecies);

  std::cout << "Reading SK-files:";
  for (iSp1 = 0; iSp1 < nSpecies; iSp1++) {
    nSK1 = angShells[iSp1].length();
    for (iSp2 = iSp1; iSp2 < nSpecies; iSp2++) {
      nSK2 = angShells[iSp2].length();
      //allocate(skData12(nSK2, nSK1)) No need allocate space, right?
      //allocate(skData21(nSK1, nSK2))
      ind = 1;
      for (iSK1 = 0; iSK1 < nSK1; iSK1++) {
        readRep = iSK1 == 0 && iSK2 == 0;
        //readAtomic = iSp1 == iSp2 && iSK1 == iSK2
        //call get(skFiles(iSp2, iSp1), fileName, ind)
      }
    }
  }

}

#endif