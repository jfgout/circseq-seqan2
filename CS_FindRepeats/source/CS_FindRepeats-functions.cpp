#include "CS_FindRepeats-functions.hpp"


std::string trimReadID(seqan::CharString &id){
  std::stringstream ss;
  std::string s;
  ss << id;
  ss >> s;

  int pos = s.find(' ');
  if( pos == std::string::npos ){
    pos = s.find('\t');
    if( pos == std::string::npos ){
      pos = s.find('/');
    }
  }

  if( pos == std::string::npos ){
    return s;
  }

  std::string st = s.substr(0, pos);
  return st;
}



int getEditDistance(std::string &s, std::string & motif, int startPos, int maxEdit){
  int i, j;
  int nbEdits = 0;
  int posMax = startPos + motif.size();
  for(i=startPos, j=0 ; i < posMax && nbEdits <= maxEdit ; i++, j++){
    if(s[i] != motif[j]){ nbEdits++; }
  }
  return nbEdits;
}


int findFirstOccurence(std::string &s, std::string &motif, int startPos, int maxEdit, int *nbEditsBest){

  int motifLg = motif.size();
  int sLg = s.size();
  int maxStart = sLg - motifLg;
  int pos;
  int editDistance;

  int bestStart = -1;
  int bestEdit = maxEdit;
  
  for(pos = startPos ; pos <= maxStart ; pos++){
    editDistance = getEditDistance(s, motif, pos, maxEdit);
    if( editDistance <= maxEdit ){
      (*nbEditsBest) = editDistance;
      return pos;
    }
  }

  (*nbEditsBest) = -1;
  return -1;
}



int getConsensus(std::string &seq, std::string &qual, int start, int csLg, std::string &consensus, std::string &consensusQual){
  std::vector<CS_ConsensusPosition> vp(csLg);
  return 0;
}



/*! \fn 
  \brief Processes a batch of seqan entries (because read/pairs are read in batches from the fastq files)  and fills a vector of CS_Repeats.
 */
void findRepeat_processEntry(std::vector<CS_Repeats> &vcs,
			     int i,
			     seqan::StringSet<seqan::CharString> &ids1,
			     seqan::StringSet<seqan::CharString> &ids2,
			     seqan::StringSet<seqan::Dna5String> &seqs1,
			     seqan::StringSet<seqan::Dna5String> &seqs2,
			     seqan::StringSet<seqan::CharString> &quals1,
			     seqan::StringSet<seqan::CharString> &quals2,
			     int minRepeatSize,
			     int maxEdits,
			     double maxDivergence,
			     int minQuality,
			     BasesIndices &bi,
			     int illuminaMinScore, int illuminaMaxScore,
			     bool paired
			     ){


  std::string sid1;
  std::string sRead, sPair, qRead, qPair;
  std::stringstream ssid;
  std::stringstream ss1;
  std::stringstream ss2;
  std::stringstream sq1;
  std::stringstream sq2;

  sid1 = trimReadID(ids1[i]);

  // Converting sequan Dna5String and CharString vairables into C++ strings
  ss1 << seqs1[i];
  ss1 >> sRead;
  if( paired == true ){
    ss2 << seqs2[i];
    ss2 >> sPair;
    sq2 << quals2[i];
    sq2 >> qPair;
  }
  
  sq1 << quals1[i];
  sq1 >> qRead;
  
  CS_Repeats csr(sid1,
		 sRead, sPair, qRead, qPair,
		 minRepeatSize, maxEdits, maxDivergence, minQuality,
		 bi, illuminaMinScore, illuminaMaxScore,
		 paired);
  //csr.getConsensus(consensusSeq, consensusQual, &avgId, illuminaMinScore, illuminaMaxScore);
  vcs[i] = csr;
}
