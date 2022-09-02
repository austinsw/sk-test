#ifndef T_SLATER_H
#define T_SLATER_H

#include <iostream>
#include <vector>
#include "t_orbitals.h"
#include "t_pair_repulsive_item.h"
#include "t_slako_cont.h"

// Slater-Koster data

class TSlater {

public:
  TSlater(int nSpecies): skHamCont(nSpecies), skOverCont(nSpecies){}
  std::vector<std::vector<double>> skSelf;
  std::vector<std::vector<double>> skHubbU;
  std::vector<std::vector<double>> skOcc;
  std::vector<double> mass;

  TSlakoCont skHamCont;
  TSlakoCont skOverCont;
  TOrbitals orb;
  std::vector<std::vector<TPairRepulsiveItem>> pairRepilsives; //pairRepulsives(:,:)
  
};

#endif