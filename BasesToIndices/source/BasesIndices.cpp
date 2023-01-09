#include "BasesIndices.hpp"



BasesIndices::BasesIndices(void){
  int i;
  for(i=0 ; i<NB_BASES ; i++){
    this->BaseToIndice_[this->tabBases_[i]] = i;
  }  
}


int BasesIndices::getIndice(char base){
  if( base<0 || base > 'z'+'Z' ){ return -1; }
  return this->BaseToIndice_[base];
}

int BasesIndices::getIndiceCpt(int i){
  switch(i){
  case 0:
    return 1;
  case 1:
    return 0;
  case 2:
    return 3;
  case 3:
    return 2;
  default:
    return -1;
  }
  return -1;
}

char BasesIndices::getBaseCpt(char base){

  switch(base){
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  default:
    return 'N';
  }
  return 'N';
}

char BasesIndices::getBaseFromIndice(int i){
  if(i<0 || i>=NB_BASES){ return UNKNOWN_BASE; }
  return this->tabBases_[i];
}

int BasesIndices::getStrandIndice(char strand){
  switch(strand){
  case '+':
    return 0;
  case '-':
    return 1;
  default:
    return -1;
  }
  return -1;
}

char BasesIndices::getStrandFromIndice(int i){

  switch(i){
  case 0:
    return '+';
  case 1:
    return '-';
  default:
    return '?';
  }
  return '?';  
}
