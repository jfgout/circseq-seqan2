#ifndef _CS_CONSENSUS_HPP_
#define _CS_CONSENSUS_HPP_

#include "BasesIndices.hpp"

#include <string>
#include <vector>
#include <map>

#define INIT_DEPTH 20 // How much size to reserve for the vector at each position in the consensus.

class CS_ConsensusPosition{
public:
  CS_ConsensusPosition();
  bool init(std::vector<char> &vCalls, std::vector<char> &vQuals);
  bool init(std::vector<char> &vCalls, std::vector<char> &vQuals, char cs, int csQual, int altQual);
  bool set(BasesIndices &BI);
  char getConsensus(void);
  int getCsQual(void);
  int getAltQual(void);
  int getDepth(void);
  bool getParameters(char *cs, int *csQual, int *altQual, int *depth);
  bool getParameters(char *cs, int *csQual, int *altQual, int *depth, int *nbDivergent);
  bool add(char call, int qual);
  std::vector<char> getCalls(void);
  std::vector<char> getQualities(void);
  bool getNbMatchAndMismatch(int *nbM, int *nbMM);
  bool print(std::ostream &os);
private:
  std::vector<char> vCalls_;
  std::vector<char> vQuals_; // Already scaled ?
  int depth_;
  char cs_; // The consensus base
  int csQual_; // The sum of quality scores for the consensus
  int altQual_; // The sum of quality scores for non-consensus bases
};


class CS_Consensus{
public:
  CS_Consensus(void);
  //CS_Consensus( std::vector<std::string> &vSeq, std::vector<std::string> &vQual);
  //CS_Consensus( std::vector<std::string> &vSeq, std::vector<std::string> &vQual, int csLg);
  CS_Consensus( std::vector<std::string> &vSeq, std::vector<std::string> &vQual, int csLg, std::vector<int> &vPosFirstBreakpoint, BasesIndices &BI, int illuminaMinScore);
  //CS_Consensus( std::vector<std::string> &vSeq, std::vector<std::string> &vQual, std::string debreakedCS);

  int getLength(void);
  int getDepth(int pos);
  bool getParameters(int pos, char *cs, int *csQual, int *altQual, int *depth);
  double getPrctMismatchesInsideConsensus(int *nbMatch, int *nbMismatch);
  
  std::string getConsensus(void);
  bool getConsensus(std::string &seq, std::string &qual);
  
  bool getPosition(std::vector<char> &vCalls, std::vector<char> &vQuals);
  bool printDetails(std::ostream &os, BasesIndices &BI);
  //bool makeMapOfConsensusPositions(std::map<int, CS_ConsensusPosition> &m, BasesIndices &BI);
  bool setMapOfConsensusPositions(BasesIndices &BI, int illuminaMinScore);  
private:
  std::vector<std::string> vSeq_;
  std::vector<std::string> vQual_;
  int nbSeq_;
  int csLg_;
  std::string consensusSeq_;
  std::vector<int> vPosFirstBreakPoint_;
  std::vector<CS_ConsensusPosition> vPos_;
  //std::map<int, CS_ConsensusPosition> mPos_;
};


//bool makeMapOfConsensusPositions(std::map<int , CS_ConsensusPosition> &m, CS_Consensus &cs);

#endif
