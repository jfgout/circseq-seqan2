#include "gff-utils.hpp"

#include <seqan/basic.h>
#include <seqan/gff_io.h>

#include <iterator>



void getScaffoldPerTranscript(std::map<std::string, std::string> &m, const char *fName){

  seqan::GffFileIn gffIn(seqan::toCString(fName));

  std::map<std::string , std::string>::iterator itm;

  int nbRecords = 0;
  
  // Copy the file record by record.
  seqan::GffRecord record;
  while (!seqan::atEnd(gffIn)){
    seqan::readRecord(record, gffIn);

    int nbTags = length(record.tagNames);
    bool done = false;
    for(int i=0 ; done == false && i<nbTags ; i++){
      std::string sTagName(seqan::toCString(record.tagNames[i]));
      if( sTagName == "transcript_id" ){
	std::string tID(seqan::toCString(record.tagValues[i]));
	itm = m.find(tID);
	if( itm == m.end() ){
	  std::string scaffoldID(toCString(record.ref));
	  m[tID] = scaffoldID;
	}
      }
    }

    //nbRecords++;
    //if( nbRecords%100000 == 0 ){
    //std::cerr<<nbRecords<<" ...\n";
    //}
    
  }
}


void getScaffoldAndStrandPerTranscript(std::map<std::string, std::pair<std::string , char> > &m, const char *fName){

  seqan::GffFileIn gffIn(seqan::toCString(fName));

  std::map<std::string , std::pair<std::string , char> >::iterator itm;

  int nbRecords = 0;
  
  // Copy the file record by record.
  seqan::GffRecord record;
  while (!seqan::atEnd(gffIn)){
    seqan::readRecord(record, gffIn);

    int nbTags = length(record.tagNames);
    bool done = false;
    for(int i=0 ; done == false && i<nbTags ; i++){
      std::string sTagName(seqan::toCString(record.tagNames[i]));
      if( sTagName == "transcript_id" ){
	std::string tID(seqan::toCString(record.tagValues[i]));
	itm = m.find(tID);
	if( itm == m.end() ){
	  char strand = record.strand;
	  std::string scaffoldID(toCString(record.ref));
	  std::pair<std::string, char> pp(scaffoldID, strand);
	  m[tID] = pp;
	}
      }
    }    
  }

}


bool loadMapOfTranscripts(const char *fName, std::map<std::string , GFFtranscriptCoordinates> &m){
  seqan::GffFileIn gffIn;
    if( !open(gffIn, fName)){
      std::cerr << "ERROR: Could not open GFF file: '"<<fName<<"'\n";
      return false;
    }
    return loadMapOfTranscripts(gffIn, m);
}


bool loadMapOfTranscripts(seqan::GffFileIn &gffIn, std::map<std::string , GFFtranscriptCoordinates> &m){

  std::map<std::string , std::vector<seqan::GffRecord> > mr;
  
  seqan::GffRecord record;
  while( !seqan::atEnd(gffIn) ){
    seqan::readRecord(record, gffIn);
    if( record.type == "exon" ){
      int nbTags = seqan::length(record.tagNames);
      bool done = false;
      for(int i=0 ; i<nbTags && !done ; i++){
	if( record.tagNames[i] == "Parent" ){
	  std::string parentID(seqan::toCString(record.tagValues[i]));
	  mr[parentID].push_back(record);
	  done = true;
	}
      }
    }
  }

  std::map<std::string , std::vector<seqan::GffRecord> >::iterator itm;
  for(itm=mr.begin() ; itm!=mr.end() ; itm++){
    std::string tID = itm->first;

    // Removing the leading "transcript:" in the parnt name:
    if( tID.size() > 11 && tID.substr(0, 11)=="transcript:" ){
      tID = tID.substr(11, tID.size()-11);
    }
    
    GFFtranscriptCoordinates gtc(itm->second);
    m[tID] = gtc;      
  }

  return true;
}


bool sortGffRecordIncreasingGenomic(seqan::GffRecord &g1, seqan::GffRecord &g2){
  return g1.beginPos < g2.beginPos;
}

bool sortGffRecordDecreasingGenomic(seqan::GffRecord &g1, seqan::GffRecord &g2){
  return g1.beginPos > g2.beginPos;
}
