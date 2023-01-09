#include "CS_Candidate.hpp"


CS_Candidate::CS_Candidate(){
  this->type_ = "?";
  this->gPos_ = -1;
  this->csPos_ = -1;
  //this->gBase_ = '?';
  //this->csBase_ = '?';
  this->sequenceFrom_ = "";
  this->sequenceTo_ = "";
}

CS_Candidate::CS_Candidate(int gPos, int csPos, char gBase, char csBase){
  this->type_ = "BASE_SUB";
  this->gPos_ = gPos;
  this->csPos_ = csPos;
  //this->gBase_ = gBase;
  //this->csBase_ = csBase;
  this->sequenceFrom_ = "";
  this->sequenceTo_ = "";
  this->sequenceFrom_.push_back(gBase);
  this->sequenceTo_.push_back(csBase);
}

CS_Candidate::CS_Candidate(int gPos, int csPos, std::string sequenceFrom, std::string sequenceTo, std::string errorType){
  this->type_ = errorType;
  this->gPos_ = gPos;
  this->csPos_ = csPos;
  this->sequenceFrom_ = sequenceFrom;
  this->sequenceTo_ = sequenceTo;
}


void CS_Candidate::print(std::string &readID, std::string &chrID, int csLg, std::ostream &os, char tStrand, bool mappedRC, int oneBased){
  os<<this->type_<<'\t'
    <<readID<<'\t'
    <<this->csPos_<<'\t'
    <<csLg<<'\t'
    <<chrID<<'\t'
    <<this->gPos_ + oneBased<<'\t'
    <<this->sequenceFrom_<<'\t'
    <<this->sequenceTo_<<'\t'
    <<tStrand<<'\t'
    <<mappedRC
    <<'\n';
}
