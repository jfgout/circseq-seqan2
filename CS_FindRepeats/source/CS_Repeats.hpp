#ifndef _CS_REPEATS_HPP_
#define _CS_REPEATS_HPP_

#include <iostream>

#include "CS_Consensus.hpp"
#include "CS_definitions.hpp"
//#include "CS_FindRepeats-functions.hpp"


#define MAX_EDITS 2

class CS_Repeats{
public:
  CS_Repeats(void);
  ~CS_Repeats(void);
  CS_Repeats(std::string id, std::string read, std::string pair, std::string readQual, std::string pairQual,
	     int minSize, int maxEdits, double maxDivergence, int minQuality,
	     BasesIndices &bi, int illuminaMinScore, int illuminaMaxScore,
	     bool paired);
  bool buildConsensus(std::string id, std::string read, std::string pair, std::string readQual, std::string pairQual,
		      int minSize, int maxEdits, double maxDivergence, int minQuality,
		      BasesIndices &bi, int illuminaMinScore, int illuminaMaxScore,
		      bool paired);
  bool addSequence(std::string seq, std::string qual, int csLg, int start, BasesIndices &bi, int illuminaMinScore, int minQuality);
  bool getConsensus(std::string &seq, std::string &qual, double *avgId, int illuminaMinScore, int illuminaMaxScore);
  int getPosStartInMate(void);
  bool hasRepeats(void);
  std::string getId(void);
  void print(std::ostream &os, double *avgDiv, int illuminaMinScore, int illuminaMaxScore);
  void printDetails(std::ostream &os, double *avgDiv, int illuminaMinScore, int illuminaMaxScore);
private:
  bool isRepeat_;
  int posStartInMate_;
  std::string id_;
  std::vector<CS_ConsensusPosition> vPos_;
};


#endif
