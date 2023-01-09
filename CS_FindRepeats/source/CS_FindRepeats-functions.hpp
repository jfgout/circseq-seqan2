#ifndef _FIND_REPEATS_FUNCTIONS_HPP_
#define _FIND_REPEATS_FUNCTIONS_HPP_

#include <iostream>
#include <seqan/seq_io.h>

#include "CS_Consensus.hpp"
#include "CS_Repeats.hpp"

std::string trimReadID(seqan::CharString &id);
int getEditDistance(std::string &s, std::string & motif, int startPos, int maxEdit);
int findFirstOccurence(std::string &s, std::string &motif, int startPos, int maxEdit, int *nbEditsBest);

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
			     bool paired);


#endif
