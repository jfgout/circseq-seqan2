#ifndef _BASES_INDICES_HPP_
#define _BASES_INDICES_HPP_

#define NB_STRANDS 2
#define NB_BASES 4

#define UNKNOWN_BASE 'N'

class BasesIndices{
public:
  BasesIndices();
  int getIndice(char base);
  char getBaseFromIndice(int i);
  int getIndiceCpt(int i);
  char getBaseCpt(char base);
  int getStrandIndice(char strand);
  char getStrandFromIndice(int i);
  //static void init(void);
  char tabBases_[NB_BASES] = {'A', 'T', 'C', 'G'};
  char tabBasesCpt_[NB_BASES] = {'T', 'A', 'G', 'C'};
  int BaseToIndice_['z'+'Z'] = {-1};
  int illuminaBase33 = 33;
};


#endif
