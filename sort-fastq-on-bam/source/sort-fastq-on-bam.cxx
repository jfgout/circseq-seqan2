#include "unistd.h"

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <map>

#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

void my_truncateID(seqan::CharString &s){
  int i;
  int lg = length(s);
  bool trimIt = false;
  for(i = 0 ; trimIt == false && i<lg ; i++){
    if( s[i]==' ' || s[i]=='\t' || s[i]=='/')
      trimIt = true;
  }
  if( trimIt==true ){
    resize(s, i-1);
  }
}


int main(int argc, char *argv[]){


  seqan::ArgumentParser parser(argv[0]);
  bool pairedFile = false;

  //addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("b", "bam", "The sam/bam file to use as input",
                                          seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("1", "r1", "The fastq file containing read sequences",
                                          seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("2", "r2", "The fastq file containing the paired read sequences",
                                          seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("o", "out1", "output file for the sorted R1 reads",
                                          seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("p", "out2", "output file for the sorted R2 reads",
                                          seqan::ArgParseArgument::STRING, "TEXT"));

  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  seqan::CharString bamFileName;
  seqan::getOptionValue(bamFileName, parser, "bam");

  seqan::CharString r1File;
  seqan::getOptionValue(r1File, parser, "r1");

  
  seqan::CharString r2File;
  if( isSet(parser, "r2") ){
    pairedFile = true;
    seqan::getOptionValue(r2File, parser, "r2");
  } else {
    r2File = "/dev/null";
  }

  seqan::CharString out1File;
  seqan::getOptionValue(out1File, parser, "out1");
  
  seqan::CharString out2File;
  if( isSet(parser, "out2") ){
    seqan::getOptionValue(out2File, parser, "out2");
  } else {
    out2File = "toDelete.fastq.gz";
  }
  
  // The fastq output will be gzip compressed on the fly with Boost:
  std::ofstream ofs1(toCString(out1File), std::ios_base::out | std::ios_base::binary);  
  std::ofstream ofs2(toCString(out2File), std::ios_base::out | std::ios_base::binary);

  boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf1;
  boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf2;
  
  outbuf1.push(boost::iostreams::gzip_compressor());
  outbuf2.push(boost::iostreams::gzip_compressor());

  outbuf1.push(ofs1);
  outbuf2.push(ofs2);

  std::ostream out1(&outbuf1);
  std::ostream out2(&outbuf2);
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // First part: read all the fastq reads and store them in mReads (where the key is the read ID for fast access)
  //  
  std::map<seqan::CharString, int> mReads;
  
  int BATCH_SIZE = 10000;
  seqan::StringSet<seqan::CharString> ids1;
  seqan::StringSet<seqan::Dna5String> seqs1;
  seqan::StringSet<seqan::CharString> quals1;

  seqan::StringSet<seqan::CharString> ids2;
  seqan::StringSet<seqan::Dna5String> seqs2;
  seqan::StringSet<seqan::CharString> quals2;

  seqan::SeqFileIn r1FileIn(seqan::toCString(r1File));
  seqan::SeqFileIn r2FileIn(seqan::toCString(r2File));

  std::cerr<<"Loading all the reads in memory (from '"<<r1File<<"' ...";
  seqan::readRecords(ids1, seqs1, quals1, r1FileIn);
  std::cerr<<" done.\n";

  if( pairedFile == true ){
    std::cerr<<"Loading all the reads in memory (from '"<<r2File<<"' ...";  
    seqan::readRecords(ids2, seqs2, quals2, r2FileIn);
    std::cerr<<" done.\n";
  }
  
  int nbRecordsRead = seqan::length(ids1);

  std::cerr<<"Indexing the reads ...";
  for(int i=0 ; i<nbRecordsRead ; i++){
    //std::cerr<<"\n'"<<ids1[i]<<"' --> '";
    my_truncateID(ids1[i]);

    if( pairedFile == true ){
      my_truncateID(ids2[i]);
    }
    //std::cerr<<ids1[i]<<"'\n";
    if( pairedFile == true &&  ids1[i] != ids2[i] ){
      std::cerr<<"ERROR, THE TWO READ FILES ARE OUT OF SYNC AT POSITION "<<i<<" ('"<<ids1[i]<<"' - '"<<ids2[i]<<"')\n";
      std::cerr<<"ABORTING THE PROGRAM.\n";
      return -1;
    }
    mReads[ids1[i]] = i;
  }
  std::cerr<<" done.\n";  

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Second part: read the bam file and output the corresponding reads in the same order
  //
  seqan::BamFileIn bamFileIn(seqan::toCString(bamFileName));

  // Need to read the header before the records...
  seqan::BamHeader header;
  seqan::readHeader(header, bamFileIn);
  
  seqan::BamAlignmentRecord record;
  std::map<seqan::CharString , int>::iterator itm;
  std::cerr<<"Sorting and outputing sorted read ...\n";
  while (!seqan::atEnd(bamFileIn)){
    seqan::readRecord(record, bamFileIn);
    itm = mReads.find(record.qName);
    if( itm == mReads.end() ){
      std::cerr<<"ERROR: '"<< record.qName <<"' NOT FOUND! ABORTING...\n";
      return -2;
    } else {
      int rPos = itm->second;
      out1<<"@"<<record.qName<<"\n"
	  <<seqs1[rPos]
	  <<"\n+\n"
	  <<quals1[rPos]
	  <<"\n";

      if( pairedFile == true ){
	out2<<"@"<<record.qName<<"\n"
	    <<seqs2[rPos]
	    <<"\n+\n"
	    <<quals2[rPos]
	    <<"\n";
      }
    }
  }

  boost::iostreams::close(outbuf1);
  boost::iostreams::close(outbuf2);
  ofs1.close();
  ofs2.close();
  
  return 0;
}

