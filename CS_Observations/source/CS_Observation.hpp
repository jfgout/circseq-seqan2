#ifndef _CS_OBSERVATION_HPP_
#define _CS_OBSERVATION_HPP_

#include <iostream>

#include "BasesIndices.hpp"


class CS_Observation{
public:
  CS_Observation();
  bool addObservation(char base, int nb, BasesIndices &BI);
  bool addObservation(int ibase, int nb);
  int getObservation(char base, BasesIndices &BI);
  int getObservation(int iBase);
  void print(std::ostream &os);
  bool isEmpty(void);
private:
  int nb_[NB_BASES];
};

CS_Observation combineObservationsFromTwoStrands(CS_Observation &obsPlus, CS_Observation &obsMinus, BasesIndices &BI);

#endif
