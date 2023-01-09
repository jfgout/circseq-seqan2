//###############################################################################################
//
// This program computes the number of BamAlignmentRecords per read


#include "unistd.h"

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include <map>
#include <vector>
#include <iterator>

int main(int argc, char const *argv[]){

  int i;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Argument parsing.
  seqan::ArgumentParser parser("refine");

  //addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption(
        "k", "kallisto-bam", "The sam/bam file produced by Kallisto",
        seqan::ArgParseArgument::STRING, "TEXT"));
  
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;
  
  seqan::CharString bamFileName;
  seqan::getOptionValue(bamFileName, parser, "kallisto-bam");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Open bam file.
  seqan::BamFileIn bamFileIn;
  if (!seqan::open(bamFileIn, seqan::toCString(bamFileName))){
    std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
    return 1;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read header.
  seqan::BamHeader header;
  seqan::readHeader(header, bamFileIn);
  typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;
  TBamContext const & bamContext = context(bamFileIn);
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Read records.
  seqan::BamAlignmentRecord record;
  int nbTotalBamRecordsRead = 0;
  bool atEndOfBamFile = false;
  std::string previousID = "";
  int nbBamReadsRead = 0;
  
  while( atEndOfBamFile == false ){

    readRecord(record, bamFileIn);
    nbTotalBamRecordsRead++;
    if( nbTotalBamRecordsRead % 1000000 == 0 ){
      fprintf(stderr, "%d bam records read ...\n", nbTotalBamRecordsRead);
    }
      
    std::string readID(seqan::toCString(record.qName));
    if( readID != previousID ){
      if( nbBamReadsRead > 0 ){
	fprintf(stdout, "%s\t%d\n", previousID.c_str(), nbBamReadsRead);
      }
      previousID = readID;
      nbBamReadsRead = 0;
    }

    nbBamReadsRead++;

    atEndOfBamFile = seqan::atEnd(bamFileIn);
  }
  fprintf(stdout, "%s\t%d\n", seqan::toCString(record.qName), nbBamReadsRead);   
    
  return 0;
}


