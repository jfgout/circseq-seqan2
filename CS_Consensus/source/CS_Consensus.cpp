#include "CS_Consensus.hpp"

#include <iostream>

CS_ConsensusPosition::CS_ConsensusPosition(){
  this->cs_ = 'N';
  this->csQual_ = 0;
  this->altQual_ = 0;
  this->depth_ = 0;
  this->vCalls_.resize(INIT_DEPTH);
  this->vQuals_.resize(INIT_DEPTH);
  //std::cerr<<"RESERVING...\n";
}


bool CS_ConsensusPosition::init(std::vector<char> &vCalls, std::vector<char> &vQuals){

  int nb = vCalls.size();
  if( nb != vQuals.size() ){  return false; }
  vCalls.swap(this->vCalls_);
  vQuals.swap(this->vQuals_);
  /*
  // Let's start by clearing the vectors in case they already had some data
  std::vector<char>().swap(this->vCalls_);
  std::vector<int>().swap(this->vQuals_);

  this->vCalls_.resize(nb);
  this->vQuals_.resize(nb);

  for(int i=0 ; i<nb ; i++){
    this->vCalls_[i] = vCalls[i];
    this->vQuals_[i] = vQuals[i];
  }
  */  
  return true;
}


bool CS_ConsensusPosition::set(BasesIndices &BI){
  int tabQuals[NB_BASES];
  for(int iBase = 0 ; iBase < NB_BASES ; iBase++){
    tabQuals[iBase] = 0;
  }
  
  //int depth = this->vCalls_.size();
  int depth = this->depth_;
  for(int iDepth = 0 ; iDepth < depth ; iDepth++){
    int iBase = BI.getIndice(this->vCalls_[iDepth]);
    tabQuals[iBase] += this->vQuals_[iDepth];
  }
  int qMax = -1;
  int whichMax = -1;
  for(int i=0 ; i<NB_BASES ; i++){
    if( tabQuals[i] > qMax ){
      qMax = tabQuals[i];
      whichMax = i;
    }
  }
  this->cs_ = BI.getBaseFromIndice(whichMax);
  this->csQual_ = qMax;
  int altQ = 0;
  for(int i=0 ; i<NB_BASES ; i++){
    if(i != whichMax){ altQ += tabQuals[i]; }
  }
  this->altQual_ = altQ;
  return true;
}


bool CS_ConsensusPosition::init(std::vector<char> &vCalls, std::vector<char> &vQuals, char cs, int csQual, int altQual){
  this->init(vCalls, vQuals);
  this->cs_ = cs;
  this->csQual_ = csQual;
  this-> altQual_ = altQual;
  return true;
}


char CS_ConsensusPosition::getConsensus(void){ return this->cs_; }

int CS_ConsensusPosition::getCsQual(void){ return this->csQual_; }

int CS_ConsensusPosition::getAltQual(void){ return this->altQual_; }

int CS_ConsensusPosition::getDepth(void){ return this->depth_; }

bool CS_ConsensusPosition::getParameters(char *cs, int *csQual, int *altQual, int *depth){
  (*cs) = this->cs_;
  (*csQual) = this->csQual_;
  (*altQual) = this->altQual_;
  //(*depth) = this->vCalls_.size();
  (*depth) = this->depth_;
  return true;
}

bool CS_ConsensusPosition::getParameters(char *cs, int *csQual, int *altQual, int *depth, int *nbDivergent){
  (*nbDivergent) = 0;
  bool res = this->getParameters(cs, csQual, altQual, depth);
  if( res == false ){ return res; }
  for(int i=0 ; i<this->depth_ ; i++){
    if( this->vCalls_[i] != (*cs) ){ (*nbDivergent)++; }
  }
  return true;
}

bool CS_ConsensusPosition::add(char call, int qual){
  //std::cerr<<"adding at position: "<<this->depth_<<" of vector of size: "<<this->vCalls_.size()<<'\n';
  //std::cerr.flush();
  this->vCalls_[this->depth_] = call;
  this->vQuals_[this->depth_] = qual;
  this->depth_ = this->depth_ + 1;
  if( this->depth_ >= this->vCalls_.size() ){
    this->vCalls_.resize(this->depth_ + INIT_DEPTH);
    this->vQuals_.resize(this->depth_ + INIT_DEPTH);
    //std::cerr<<"WARNING: DEPTH TOO DEEP!\n";
  }
  /*
  this->vCalls_.push_back(call);
  this->vQuals_.push_back(qual);
  return true;
  */
  return true;
}


bool CS_ConsensusPosition::print(std::ostream &os){
  int nb = this->vCalls_.size();
  for(int i=0 ; i<nb ; i++){
    os<<this->vCalls_[i];
  }
  os<<"\t("
      <<this->csQual_
      <<" - "
      <<this->altQual_
      <<")";

  return true;
}

std::vector<char> CS_ConsensusPosition::getCalls(void){
  return this->vCalls_;
}

std::vector<char> CS_ConsensusPosition::getQualities(void){
  return this->vQuals_;
}

bool CS_ConsensusPosition::getNbMatchAndMismatch(int *nbM, int *nbMM){
  (*nbM) = 0;
  (*nbMM) = 0;
  for(int i=0 ; i<this->depth_ ; i++){
    if( this->vCalls_[i] == this->cs_ ){
      (*nbM)++;
    } else {
      (*nbMM)++;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


CS_Consensus::CS_Consensus(std::vector<std::string> &vSeq, std::vector<std::string> &vQual, int csLg, std::vector<int> &vPosFirstBreakPoint, BasesIndices &BI, int illuminaMinScore){
  this->vSeq_ = vSeq;
  this->vQual_ = vQual;
  this->vPosFirstBreakPoint_ = vPosFirstBreakPoint;
  this->csLg_ = csLg;
  this->setMapOfConsensusPositions(BI, illuminaMinScore);
}

std::string CS_Consensus::getConsensus(void){
  return this->consensusSeq_;
}

int CS_Consensus::getLength(void){
  return this->csLg_;
}

int CS_Consensus::getDepth(int pos){
  if( pos >= this->csLg_ ){ return -1; }
  return this->vPos_[pos].getDepth();
}

bool CS_Consensus::getParameters(int pos, char *cs, int *csQual, int *altQual, int *depth){
  if( pos<0 || pos >= this->csLg_ ){ return -1; }
  return this->vPos_[pos].getParameters(cs, csQual, altQual, depth);
}

bool CS_Consensus::setMapOfConsensusPositions(BasesIndices &BI, int illuminaMinScore){

  int csLg = this->csLg_;
  int nbSeq = this->vSeq_.size();

  //CS_ConsensusPosition csp;
  this->vPos_.resize(csLg);

  for(int iSeq=0 ; iSeq < nbSeq ; iSeq++){
    int seqLength = this->vSeq_[iSeq].size();
    int posFirstBreakPoint = this->vPosFirstBreakPoint_[iSeq];
    for(int sPos=0 ; sPos<seqLength ; sPos++){
    // Convert position in sequence into position in (debreaked) consensus
      int csPos = ((sPos - posFirstBreakPoint)+2*csLg)%csLg; // Here I have to add 2*csLg because c++ modulo does not work the same as R modulo --> can return negative value!
      this->vPos_[csPos].add(this->vSeq_[iSeq][sPos], this->vQual_[iSeq][sPos] - illuminaMinScore);
    }
  }

  std::string consensusSeq(csLg, '?');
  for(int i=0 ; i<csLg ; i++){
    this->vPos_[i].set(BI);
    consensusSeq[i] = this->vPos_[i].getConsensus();
  }
  this->consensusSeq_ = consensusSeq;

  // Alternative (might be faster?): for each position in the consensus, find all corresponding positions in the sequences and construct the corresponding CS_ConsensusPosition
  return true;
}


bool CS_Consensus::getConsensus(std::string &seq, std::string &qual){
  if( this->csLg_ == 0 ){
    return false;
  }

  seq = this->consensusSeq_;
  qual = seq;
  int lg = seq.size();
  for(int i=0 ; i<lg ; i++){
    int iQual = this->vPos_[i].getCsQual();
    if( iQual > 42 ){ iQual = 42; }
    char cQual = (char)(iQual + 33);
    qual[i] = cQual;
  }
  return true;
}

bool CS_Consensus::printDetails(std::ostream &os, BasesIndices &BI){

  for(int i=0 ; i<this->csLg_ ; i++){
    os<<i<<'\t';
    this->vPos_[i].print(os);
    os<<'\n';
  }
  
  /*
  std::map<int, CS_ConsensusPosition>::iterator itm;
  for(itm=this->mPos_.begin() ; itm!=this->mPos_.end() ; itm++){
    os<<itm->first<<'\t';
    (itm->second).print(os);
    os<<'\n';
  }
  */
  return true;  
}



double CS_Consensus::getPrctMismatchesInsideConsensus(int *nbMatch, int *nbMismatch){
  (*nbMatch) = 0;
  (*nbMismatch) = 0;
  for(int i=0 ; i<this->csLg_ ; i++){
    int nbM, nbMM;
    this->vPos_[i].getNbMatchAndMismatch(&nbM, &nbMM);
    (*nbMatch) += nbM;
    (*nbMismatch) += nbMM;
  }
  int total = (*nbMatch) + (*nbMismatch);
  double res = (double)((double)(*nbMismatch)/(double)total);
  return res;
}
