#ifndef _CS_CANDIDATE_HPP_
#define _CS_CANDIDATE_HPP_

#include <iostream>

class CS_Candidate{
public:
  CS_Candidate();
  CS_Candidate(int gPos, int csPos, char gBase, char csBase); // Use this for base subs
  CS_Candidate(int gPos, int csPos, std::string sequenceFrom, std::string sequenceTo, std::string errorType); // Use this for indels
  void print(std::string &readID, std::string &chrID, int csLg, std::ostream &os, char tStrand, bool mappedRC, int oneBased);
private:
  int gPos_;
  int csPos_;
  //char gBase_;
  //char csBase_;
  std::string type_;
  std::string sequenceFrom_;
  std::string sequenceTo_;
};


#endif
