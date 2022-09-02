#ifndef T_SLAKO_CONT_H
#define T_SLAKO_CONT_H

#include <iostream>
#include <vector>

class TSlakoEqGrid{};

class TSlaKo_ {
  int iType = 0;
  TSlakoEqGrid pSlakoEqGrid;
};


class TSlakoCont {

  std::vector<std::vector<TSlaKo_>> slakos; //slakos(:,:)
  int nSpecies;
  int mInt;
  double cutoff;
  bool tDataOK;
  bool tInit = false;

// SlakoCont_init
public:
  TSlakoCont(int nSpecies): nSpecies(nSpecies), mInt(0) {
    this->slakos.resize(nSpecies, std::vector<TSlaKo_>(nSpecies));
    //this->mInt = 0;
    this->cutoff = 0.0;
    this->tDataOK = false;
    this->tInit = true;
  }

};

#endif