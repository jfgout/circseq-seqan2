#include "CS_FindRepeats-functions.hpp"
#include "CS_Repeats.hpp"

CS_Repeats::~CS_Repeats(void){
  this->vPos_.clear();
}


CS_Repeats::CS_Repeats(void){
  this->isRepeat_ = false;
  this->posStartInMate_ = -1;
  this->id_ = "";
}

CS_Repeats::CS_Repeats(std::string id,
		       std::string read, std::string pair, std::string readQual, std::string pairQual,
		       int minSize, int maxEdits, double maxDivergence, int minQuality,
		       BasesIndices &bi, int illuminaMinScore, int illuminaMaxScore,
		       bool paired){
  this->buildConsensus(id, read, pair, readQual, pairQual, minSize, maxEdits, maxDivergence, minQuality, bi, illuminaMinScore, illuminaMaxScore, paired);
}


bool CS_Repeats::buildConsensus(std::string id,
				std::string read, std::string pair, std::string readQual, std::string pairQual,
				int minSize, int maxEdits, double maxDivergence, int minQuality,
				BasesIndices &bi, int illuminaMinScore, int illuminaMaxScore,
				bool paired){

  this->isRepeat_ = false;
  this->id_ = id;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Let's see if the read has at least one repeat:
  std::string stmp = read.substr(0, minSize); // Extracting the seed
  int nbEdits = -1;
  int pos = findFirstOccurence(read, stmp, minSize, maxEdits, &nbEdits); // Searching for a first repeat of the seed in the read
  // Possible improvement: allow for partial repeat. For example, read of size 200 with repeat of size 120 (partial repeat will be 80nt long)
  // Need to check it that's already the case.
  if( pos <= 0 )
    return false;

  int csLg = pos;
  this->vPos_.resize(csLg);
  this->addSequence(read, readQual, csLg, 0, bi, illuminaMinScore, minQuality);

  // In case I'm looking at single-end data (pair and pairQual are empty):
  if( paired == false ){
    this->isRepeat_ = true;
    return true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////
  // Searching for consensus in the pair

  // I start by obtaining a clean consensus from the read
  std::string tmpCs, tmpQual;
  double avgDiv;
  this->getConsensus(tmpCs, tmpQual, &avgDiv, illuminaMinScore, illuminaMaxScore);
  if( avgDiv > maxDivergence ){ return false; }

  // Looking for the first occurence of the consensus in the mate:
  pos = findFirstOccurence(pair, tmpCs, 0, maxEdits, &nbEdits);
  //std::cerr<<"Searching for: '"<<tmpCs<<"' in '"<<pair<<"'\n";
  //std::cerr<<"FirstOccurence in mate: "<<pos<<std::endl;
  if( pos < 0 ){ return false; }
  if( pos > csLg ){
    pos = pos % csLg;
  }
  this->posStartInMate_ = pos;
  this->addSequence(pair, pairQual, csLg, pos, bi, illuminaMinScore, minQuality);
  
  this->isRepeat_ = true;
  return true;
}


bool CS_Repeats::addSequence(std::string seq, std::string qual, int csLg, int start, BasesIndices &bi, int illuminaMinScore, int minQuality){
  int i, csPos;
  int seqLg = seq.size();
  for(i=0 ; i<seqLg ; i++){
    if( i<start ){
      csPos = (csLg - start) + i;
    } else {
      csPos = (i - start) % csLg;
    }
    // Only base calls that have a quality score >= to 'minQuality' (typically set to 30) are incorporated in the consensus.
    int qualityScore = qual[i] - illuminaMinScore;
    if( qualityScore >= minQuality ){
      this->vPos_[csPos].add(seq[i], qualityScore);
    }
  }

  for(i=0 ; i<csLg ; i++){
    this->vPos_[i].set(bi);
  }
  
  return true;
}


/**
 This function goes through the vector of base-calls and qualities at each position of the repeat and calls the consensus + computes the average divergence from it.
*/
bool CS_Repeats::getConsensus(std::string &seq, std::string &qual, double *avgId, int illuminaMinScore, int illuminaMaxScore){
  int csLg = this->vPos_.size();
  seq.resize(csLg, 'N');
  qual.resize(csLg, (char)illuminaMinScore);

  char cs;
  int csQual, altQual, depth, nbDivergent;
  int totalBases = 0;
  int totalDivergent = 0;
  for(int i=0 ; i<csLg ; i++){
    bool res = this->vPos_[i].getParameters(&cs, &csQual, &altQual, &depth, &nbDivergent);
    if( res == false ){ return res; }
    totalBases += depth;
    totalDivergent += nbDivergent;
    seq[i] = cs;
    int localQuality = csQual + illuminaMinScore;
    if( localQuality > illuminaMaxScore ){ localQuality = illuminaMaxScore; }
    if( localQuality < illuminaMinScore ){ localQuality = illuminaMinScore; }
    qual[i] = (char)localQuality;
  }
  if( totalBases < 1 ){ return false; }
  (*avgId) = (double)totalDivergent / (double)totalBases;
  return true;
}

bool CS_Repeats::hasRepeats(void){
  return this->isRepeat_;
}

int CS_Repeats::getPosStartInMate(void){
  return this->posStartInMate_;
}

std::string CS_Repeats::getId(void){
  return this->id_;
}

void CS_Repeats::print(std::ostream &os, double *avgDiv, int illuminaMinScore, int illuminaMaxScore){
  std::string csSeq, csQual;
  this->getConsensus(csSeq, csQual, avgDiv, illuminaMinScore, illuminaMaxScore);
  os<<this->id_<<'\n'
    <<csSeq<<'\n'
    <<"+\n"
    <<csQual
    <<'\n';
}


void CS_Repeats::printDetails(std::ostream &os, double *avgDiv, int illuminaMinScore, int illuminaMaxScore){

  std::string csSeq, csQual;
  this->getConsensus(csSeq, csQual, avgDiv, illuminaMinScore, illuminaMaxScore);
  os<<this->id_<<" : "<<(*avgDiv)
    <<'\n'
    <<csSeq
    <<'\n';
  
  int csLg = csSeq.size();
  int maxDepth = 0;
  std::map<int, std::pair< std::vector<char> , std::vector<char> > >mPos;
  for(int i=0; i<csLg ; i++){
    int depth = this->vPos_[i].getDepth();
    std::vector<char> vCalls = this->vPos_[i].getCalls();
    std::vector<char> vQuals = this->vPos_[i].getQualities();
    vCalls.resize(depth);
    vQuals.resize(depth);
    std::pair< std::vector<char> , std::vector<char> > pp(vCalls, vQuals);
    mPos[i] = pp;
    if( vCalls.size() > maxDepth ){ maxDepth = vCalls.size(); }
  }

  os<<"MaxDepth: "<<maxDepth<<'\n';
  os<<mPos[0].first.size()<<'\n';
  
  for(int depth = 0 ; depth < maxDepth ; depth++){
    for(int i=0 ; i<csLg ; i++){
      if( mPos[i].first.size() <= depth ){
	os<<' ';
      } else {
	os<<mPos[i].first[depth];
      }
    }
    os<<'\n';
  }
  
}
