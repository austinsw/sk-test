// Not finished. But seems like the original code was not doing anything either.

#ifndef T_PAIR_REPULSIVE_ITEM_H
#define T_PAIR_REPULSIVE_ITEM_H

class TPairRepulsive {

public:
  double TPairRepulsive_getCutoff(TPairRepulsive self);
  void TPairRepulsive_getValue();

};

// Wrapper for a TPairRepulsive class item
struct TPairRepulsiveItem{

  TPairRepulsive item;

};

#endif