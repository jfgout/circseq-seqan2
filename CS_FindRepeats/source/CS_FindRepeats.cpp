/*
This program finds repeats in (pairs of) reads.
It outputs, for each (pair of) read with repeats, three files:

1) A "positions" file with 3 columns: readID, size of the repeats, start location of the first full-length repeat in the mate.

2) A fastq file with the consensus derived from the repeats.

3) A log file with information about the number of reads for which the program has found a repeating pattern.

As an option, it can also output a subset (or all) the reads for which the program did not find any repeat. These can beuseful for debugging.

 */

#include <iostream>
#include <string>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <omp.h>
#include <chrono>

#include "CS_FindRepeats-functions.hpp"
#include "CS_Repeats.hpp"

using namespace std;
using namespace std::chrono;

//#define NB_READS_PER_PACKET 2000000
//#define NB_READS_PER_TEST_REVERSE 10000

int main(int argc, char *argv[]){

  seqan::ArgumentParser parser(argv[0]);

  //addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("1", "R1", "The fastq file containing read sequences", seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("2", "R2", "The fastq file containing pairs sequences", seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("s", "min-size", "The minimum repeat size (default = 25)", seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption("e", "max-edits", "The maximum number of edits in the first repeat (default = 2)", seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption("d", "max-divergence", "The maximum amount of sequence divergence between the repeats (default = 0.05)", seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
  addOption(parser, seqan::ArgParseOption("q", "min-quality", "The minimum quality score for base-calls to be included in the divergence calculation (default = 30)", seqan::ArgParseArgument::INTEGER, "INT"));
  
  addOption(parser, seqan::ArgParseOption("b", "base-name", "Base name for the output files", seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("r", "reverse", "Should the R2 sequences be reverse-complemented? T = True, F = False, G = Guess (default)", seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("t", "threads", "How many threads to use? (default = 1)", seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption("p", "packet-size", "How many reads to read per batch (default = 200,000)", seqan::ArgParseArgument::INTEGER, "INT"));

  addOption(parser, seqan::ArgParseOption("z", "gz", "Use gzip compression for the output."));

  addOption(parser, seqan::ArgParseOption("m", "min-score", "The lowest possible base-call quality score (default = 33 (ASCII: '!' - illumina 'recent' data)", seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption("x", "max-score", "The highest possible base-call quality score (default = 75 (ASCII: 'K'- illumina 'recent' data)", seqan::ArgParseArgument::INTEGER, "INT"));
  
  
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  if ( isSet(parser, "R1") == false ){
    std::cerr<<"ERROR: NO INPUT FILE PROVIDED.\n";    
    return -1;
  }

  if ( isSet(parser, "base-name") == false ){
    std::cerr<<"ERROR: No basename provided.\n";  
    return -1;
  }
  
  seqan::CharString r1FileName, r2FileName;
  bool paired = false;
  bool reverseMate = false;
  seqan::CharString cReverse("G");
  
  seqan::getOptionValue(r1FileName, parser, "R1");
  if ( isSet(parser, "R2") == true ){
    seqan::getOptionValue(r2FileName, parser, "R2");
    paired = true;
    if( isSet(parser, "reverse") == true ){
      seqan::getOptionValue(cReverse, parser, "reverse");
    }
    if( cReverse == "T" ){ reverseMate = true; }
  }


  bool gzOut = false;
  if( isSet(parser, "gz") == true ){
    gzOut = true;
  }

  // Handling the output file that will contain the readID with the size of the repeat (+ position of first occurence of the repeat in paired read if data is paired-end)
  seqan::CharString cs_baseName;
  seqan::getOptionValue(cs_baseName, parser, "base-name");
  std::string fNameResults(seqan::toCString(cs_baseName));
  fNameResults = fNameResults + "-CS.pos";
  if( gzOut == true ){
    fNameResults = fNameResults + ".gz";
  }
  std::ofstream ofsResults(fNameResults, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::output> bosResults;
  if( gzOut == true ){
    bosResults.push(boost::iostreams::gzip_compressor());
  }
  bosResults.push(ofsResults);
  std::ostream osResults(&bosResults);
  // At this point, "osResults" can be use as a stream to output the results (and it will be in gz-compressed format if this option has been specified).

  // Now the file that will contain the double consensus (BASE_NAME-DCS.fastq)
  std::string fDCS(seqan::toCString(cs_baseName));  
  fDCS = fDCS + "-DCS.fastq";
  if( gzOut == true ){ fDCS = fDCS + ".gz"; }
  std::ofstream ofsDCS(fDCS, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::output> bosDCS;
  if( gzOut == true ){ bosDCS.push(boost::iostreams::gzip_compressor()); }
  bosDCS.push(ofsDCS);
  std::ostream osDCS(&bosDCS);
  // Now, the DCS.fastq file can be written in with osDCS
  

  // The size of the "seed" to use to search for a repeat
  int minRepeatSize = 25;
  if ( isSet(parser, "min-size") == true ){
    seqan::getOptionValue(minRepeatSize, parser, "min-size");
  }

  // The maximum number of mismatches between the seed and its first match
  int maxEdits = 2;
  if( isSet(parser, "max-edits") == true ){
    seqan::getOptionValue(maxEdits, parser, "max-edits");
  }

  // Report only cases where the total divergence between all repeats is < maxDivergence
  double maxDivergence = 0.05;
  if( isSet(parser, "max-divergence") == true ){
    seqan::getOptionValue(maxDivergence, parser, "max-divergence");
  }

  // Use only base calls with quality score > minQuality in the calculation of total divergence
  int minQuality = 30;
  if( isSet(parser, "min-quality") == true ){
    seqan::getOptionValue(minQuality, parser, "min-quality");
  }
  
  int nbThreads = 1;
  if( isSet(parser, "threads") == true ){
    seqan::getOptionValue(nbThreads, parser, "threads");
  }

  // Reads are loaded in memory into packets (efficient file access) of size NB_READS_PER_PACKET
  int NB_READS_PER_PACKET = 200000;
  if( isSet(parser, "packet-size") == true ){
    seqan::getOptionValue(NB_READS_PER_PACKET, parser, "packet-size");
  }


  // These options are only in the unlikely event that the user wants to process data generated with the old illumina scoring scheme
  // (In the new scheme, scores range from '!' = 0 to 'K' = 42 while the old one ranges from '@' = 0 to 'j' = 42)
  int illuminaMinScore = (int)'!';
  if( isSet(parser, "min-score") == true ){
    seqan::getOptionValue(illuminaMinScore, parser, "min-score");
  }
  int illuminaMaxScore = (int)'K';
  if( isSet(parser, "max-score") == true ){
    seqan::getOptionValue(illuminaMaxScore, parser, "max-score");
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // END OF ARGUMENTS PARSING
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // These variables will store the packets of reads
  seqan::StringSet<seqan::CharString> ids1, ids2;
  seqan::StringSet<seqan::Dna5String> seqs1, seqs2;
  seqan::StringSet<seqan::CharString> quals1, quals2;

  // Opening the reads files
  seqan::SeqFileIn seqFileIn1, seqFileIn2;
  if (!open(seqFileIn1, seqan::toCString(r1FileName))) {
    std::cerr << "ERROR: Could not open the file '"<<r1FileName<<"'.\n";
    return 1;
  }
  if (paired == true && !open(seqFileIn2, seqan::toCString(r2FileName))) {
    std::cerr << "ERROR: Could not open the file '"<<r2FileName<<"'.\n";
    return 1;
  }


  // I read a first batch of reads (and mates if the data is paired-end)
  seqan::readRecords(ids1, seqs1, quals1, seqFileIn1, NB_READS_PER_PACKET);
  if( paired == true ){
    seqan::readRecords(ids2, seqs2, quals2, seqFileIn2, NB_READS_PER_PACKET);
    if( reverseMate == true ){
      std::cerr<<"Reversing the mate\n";
      seqan::reverseComplement(seqs2);
      seqan::reverse(quals2);
    }
  }

  stringstream ss;
  int nbRead = seqan::length(ids1); // <- This will store the number of reads that have been read. When I reach the end of the fastq file, it will be zero.
  int nbReadsTotal = nbRead;
  int nbReadsWithRepeat = 0;
  int i;
  //int illuminaMinScore = 33;
  //int illuminaMaxScore = 75;
  double avgDiv;

  BasesIndices bi;
  std::string consensusSeq, consensusQual;
  double avgId;
  std::vector<CS_Repeats> vcs;

  std::cerr<<"nbRead: "<<nbRead<<"\n";
  vcs.resize(nbRead);

  // This is the main loop which will process a batch of reads and read the next batch from the fastq file:
  int nbBatchDone = 0;
  while( nbRead > 0 ){

    if( paired == true && cReverse == "G" && nbBatchDone == 0 ){
    // If this is the first batch of reads and the user asked the program to guess whether mates shold be reverse-complemented or not.
    // In this situation, I will process the first batch with both reverse-complementation and no-reverse-complementation and see which one works best.
      
      // Trying first without reverse-complementation:
      int nbWithRepeatForward = 0;
      #pragma omp parallel for num_threads(nbThreads)
      for(i=0 ; i<nbRead ; i++){
	findRepeat_processEntry(vcs, i,
				ids1, ids2,
				seqs1, seqs2,
				quals1, quals2,
				minRepeatSize, maxEdits, maxDivergence, minQuality,
				bi, illuminaMinScore, illuminaMaxScore,
				paired);
	//vcs[i].printDetails(cout, &avgDiv, illuminaMinScore, illuminaMaxScore);
	if( vcs[i].hasRepeats() == true ){
	  vcs[i].getConsensus(consensusSeq, consensusQual, &avgDiv, illuminaMinScore, illuminaMaxScore);
	  if( avgDiv <= maxDivergence ){
	    nbWithRepeatForward++;
	  }	  
	}
      }

      // Now, let's try with the reverse mate
      std::cerr<<"Looking at reversed mates...\n";
      seqan::StringSet<seqan::Dna5String> seqs2r;
      seqan::StringSet<seqan::CharString> quals2r;
      seqs2r = seqs2;
      quals2r = quals2;
      seqan::reverseComplement(seqs2r);
      seqan::reverse(quals2r);
      int nbWithRepeatReverse = 0;
      #pragma omp parallel for num_threads(nbThreads)
      for(i=0 ; i<nbRead ; i++){
	findRepeat_processEntry(vcs, i,
				ids1, ids2,
				seqs1, seqs2r,
				quals1, quals2r,
				minRepeatSize, maxEdits, maxDivergence, minQuality,
				bi, illuminaMinScore, illuminaMaxScore,
				paired);

	//vcs[i].printDetails(cout, &avgDiv, illuminaMinScore, illuminaMaxScore); // <- This is only for debugging!
	if( vcs[i].hasRepeats() == true ){
	  vcs[i].getConsensus(consensusSeq, consensusQual, &avgDiv, illuminaMinScore, illuminaMaxScore);	  
	  if( avgDiv <= maxDivergence ){
	    nbWithRepeatReverse++;
	  }
	} else {
	  ;
	  // Could insert some debugging code here to report the identity % of reads that failled.
	  //cerr<<ids1[i]<<endl;
	}
      }

      // Now, compare the two results:
      std::cerr<<nbWithRepeatForward<<" vs "<<nbWithRepeatReverse<<endl;
      if( nbWithRepeatForward > nbWithRepeatReverse ){
	reverseMate = false;
	// I need to reverse (complement) again to get back to the original for this first batch:
      } else {
	reverseMate = true;
	//seqs2 = seqs2r;
	//quals2 = quals2r;
      }
    } // END of the block of code that tries to guess the orientation of the paired reads.
    
    if( reverseMate == true ){
      seqan::reverseComplement(seqs2);
      seqan::reverse(quals2);
    }
    
    // Process batch (this part to bedone in parallel):
    #pragma omp parallel for num_threads(nbThreads)
    for(i=0 ; i<nbRead ; i++){
      findRepeat_processEntry(vcs, i,
			      ids1, ids2,
			      seqs1, seqs2,
			      quals1, quals2,
			      minRepeatSize, maxEdits, maxDivergence, minQuality,
			      bi, illuminaMinScore, illuminaMaxScore,
			      paired);
    }
    // END OF THE PARALLEL CODE!!!
    
    std::vector<CS_Repeats>::iterator itcs;
    for(itcs = vcs.begin() ; itcs != vcs.end() ; itcs++){
      bool isOk = itcs->hasRepeats();
      if( isOk == true ){
	itcs->getConsensus(consensusSeq, consensusQual, &avgDiv, illuminaMinScore, illuminaMaxScore);
	if( avgDiv <= maxDivergence ){
	  nbReadsWithRepeat++;
	  int csLg = consensusSeq.size();
	  string readID = itcs->getId();

	  //itcs->printDetails(cout, &avgId, illuminaMinScore, illuminaMaxScore);
	  osResults<<readID<<'\t'
		   <<csLg;
	  if( paired == true ){
	    int posStartInMate = itcs->getPosStartInMate();
	    osResults<<'\t'<<posStartInMate;
	  }
	  osResults<<'\n';

	// Outputing the DCS.fastq file:
	// (TO DO)
	  osDCS<<"@"<<readID<<'\n'
	       <<consensusSeq<<consensusSeq<<'\n'
	       <<"+\n"
	       <<consensusQual<<consensusQual
	       <<'\n';
	  
	}	
      } else {
	// Output some information in verbose mode?
	;
	//cerr<<itcs->getId()<<endl;
      }

    }

        
    std::cerr<<"Reading a new batch... (current tally: "<<nbReadsTotal<<")\n";
    //std::cerr<<"'"<<ids1[1]<<"'\t'"<<ids2[1]<<"'\n";
    //sid1 = trimReadID(ids1[1]);
    //sid2 = trimReadID(ids2[1]);
    //std::cerr<<"'"<<sid1<<"'\t'"<<sid2<<"'\n";    
    //return 0;

    // Now reading the next batch:
    seqan::clear(ids1);
    seqan::clear(seqs1);
    seqan::clear(quals1);    
    seqan::readRecords(ids1, seqs1, quals1, seqFileIn1, NB_READS_PER_PACKET);
    if( paired == true ){
      seqan::clear(ids2);
      seqan::clear(seqs2);
      seqan::clear(quals2);
      seqan::readRecords(ids2, seqs2, quals2, seqFileIn2, NB_READS_PER_PACKET);
    }
    nbRead = seqan::length(ids1);
    nbReadsTotal += nbRead;
    nbBatchDone++;
  }

  //outResults.close();
  std::cerr<<nbReadsTotal<<"\tTOTAL_READS\n";
  std::cerr<<nbReadsWithRepeat<<"\tREADS_WITH_REPEAT\n";
  
  return 0;
}



