#include "CS_Observation.hpp"



CS_Observation::CS_Observation(void){
  for(int i=0 ; i<NB_BASES ; i++){
    this->nb_[i] = 0;
  }
}

bool CS_Observation::addObservation(char base, int nb, BasesIndices &BI){
  int iBase = BI.getIndice(base);
  return this->addObservation(iBase, nb);
}


bool CS_Observation::addObservation(int iBase, int nb){
  if( iBase >= 0 && iBase < NB_BASES ){
    this->nb_[iBase] += nb;
    return true;
  }
  return false;
}


int CS_Observation::getObservation(char base, BasesIndices &BI){
  int iBase = BI.getIndice(base);
  return this->getObservation(iBase);
}

int CS_Observation::getObservation(int iBase){
  return this->nb_[iBase];
}


void CS_Observation::print(std::ostream &os){
  for(int i=0 ; i<NB_BASES ; i++){
    if(i>0){ os<<'\t'; }
    os<<this->nb_[i];
  }
}

bool CS_Observation::isEmpty(void){
  for(int i=0 ; i<NB_BASES ; i++){
    if( this->nb_[i] > 0 ){ return false; }
  }
  return true;
}

CS_Observation combineObservationsFromTwoStrands(CS_Observation &obsPlus, CS_Observation &obsMinus, BasesIndices &BI){
  CS_Observation obs = obsPlus;
  int icpt, nb;
  for(int i=0 ; i<NB_BASES ; i++){
    icpt = BI.getIndiceCpt(i);
    nb = obsMinus.getObservation(icpt);
    obs.addObservation(i, nb);
  }
  return obs;
}
