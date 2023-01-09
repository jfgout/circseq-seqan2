#ifndef _GFF_UTILS_HPP_
#define _GFF_UTILS_HPP_

// A collection of functions for deaing with GFF/GTF files, taking advantage of the seqan library.


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include <seqan/basic.h>
#include <seqan/gff_io.h>

#include "GFFtranscript.hpp"

void getScaffoldPerTranscript(std::map<std::string, std::string> &m, const char *fName);

void getScaffoldAndStrandPerTranscript(std::map<std::string, std::pair<std::string , char> > &m, const char *fName);



bool loadMapOfTranscripts(const char *fName, std::map<std::string , GFFtranscriptCoordinates> &m);
bool loadMapOfTranscripts(seqan::GffFileIn &fin, std::map<std::string , GFFtranscriptCoordinates> &m);


bool sortGffRecordIncreasingGenomic(seqan::GffRecord &g1, seqan::GffRecord &g2);
bool sortGffRecordDecreasingGenomic(seqan::GffRecord &g1, seqan::GffRecord &g2);


#endif
