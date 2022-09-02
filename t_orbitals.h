#ifndef T_ORBITALS_H
#define T_ORBITALS_H

#include <iostream>
#include <vector>
using namespace std;

// A temporary replacement for "dftbp_type_commontypes, only : TOrbitals"
// Contains information about the orbitals of the species/atoms in the system

class TOrbitals {

public:

  // Nr. of shells for each atomic species (nSpecies)
  vector<int> nShell;

  // Nr. of orbitals for each atomic species (nSpecies)
  vector<int> nOrbSpecies;

  // Nr. of orbitals for each atom (nAtom)
  vector<int> nOrbAtom;

  // Ang. momentum of the a particular l-shell on a particular species (maxval(nShell), nSpecies)
  vector<vector<int>> angShell;

  // The shell which contains the given orbital on an atom
  // (maxval(nOrbSpecies), nSpecies)
  vector<vector<int>> iShellOrb;

  // Starting pos. within the atomic block of the each of the shells of each species
  // (maxval(nShell)+1, nSpecies)
  vector<vector<int>> posShell;

  // Max. nr. of shells for any species
  int mShell;

  // Max. nr. of orbitals for any species
  int mOrb;

  // Total number of orbitals in system.
  int nOrb;
  
};

#endif