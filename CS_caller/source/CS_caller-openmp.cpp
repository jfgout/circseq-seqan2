/*
call: the program that reads the sorted/refined BAM file and makes the calls regarding observations and candidates

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

#include "gff-utils.hpp"
#include "GFFtranscript.hpp"
#include "CS_ScaffoldObservations.hpp"
#include "CS_Consensus.hpp"
#include "CS_Candidate.hpp"

#define MAX_POS_DIFF_BEFORE_WARNING 3 // If a consensus has more than this number of differences with the reference genome, I print a warning message (this could indicate a bad alignment)
//#define NB_STRANDS 2

using namespace std;
using namespace std::chrono;

int callThisConsensus(CS_Consensus &cs,
		      int nbExcludeStart, int nbExcludeEnd,
		      seqan::BamAlignmentRecord &record,
		      std::string &chrID, char tStrand, bool mappedRC, bool strandSpecific,
		      seqan::CharString &ChrSeq,
		      GFFtranscriptCoordinates &gt,
		      std::vector<int> &vTtoG,
		      CS_ScaffoldObservations &cobCalls,
		      CS_ScaffoldObservations &cobObs,
		      int minDepth,
		      int allowAlternativeCalls,
		      int maxAltQual,
		      int minConfidenceQual,
		      int nbMaxMMperCs,
		      BasesIndices &BI, int oneBased);


void loadAllPositions(char *fName, std::map<std::string, std::pair<int,int> > &mPos);

void loadAllPositions(char *fName, std::map<std::string, std::pair<int,int> > &mPos){
  std::ifstream infile(fName);
  std::string line, tID, stmp1, stmp2;
  int pos1, pos2;
  while(getline(infile, line)){
    std::istringstream iss(line);
    std::getline(iss, tID, '\t');
    std::getline(iss, stmp1, '\t');
    std::getline(iss, stmp2, '\t');
    pos1 = atoi(stmp1.c_str());
    pos2 = atoi(stmp2.c_str());
    std::pair<int , int> pp(pos1, pos2);
    mPos[tID] = pp;
  }

}


int main(int argc, char *argv[]){

  seqan::ArgumentParser parser(argv[0]);
  bool pairedFile = false;

  //addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("b", "bam", "The sam/bam file to use as input", seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("p", "pos", "The file that contains the consensus length and position in pair (generated in the first step of the pipeline)", seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("1", "R1", "The fastq file containing read sequences", seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("2", "R2", "The fastq file containing pairs sequences", seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("r", "reference-genome", "The reference genome (uncompressed fasta, pre-indexed with the corresponding .fai file)", seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption("g", "GFF", "A GFF/GTF annotation file containing the exon coordinates for all the transcripts", seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("c", "calls", "Filename where to output the list of calls made per genomic position.", seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("o", "observations", "Filename where to output the list of all observations (including low-quality ones - used to find het/polymoprhic sites)", seqan::ArgParseArgument::STRING, "TEXT"));

  
  addOption(parser, seqan::ArgParseOption("d", "min-depth", "Minimum depth in the consensus to make a call (recommended/default value = 3).", seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption("q", "min-consensus-qual", "Minimum quality of the consensus (defined as [sum of quality scores for the consensus base]-[sum of quality score for alternative bases] (recommended value: 100).", seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption("exclude-start", "exclude-start", "How many nucleotides to skip at the beginning of each consensus?.", seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption("exclude-end", "exclude-end", "How many nucleotides to skip at the end of each consensus?.", seqan::ArgParseArgument::INTEGER, "INT"));

  addOption(parser, seqan::ArgParseOption("a", "allow-alternative-calls", "Should positions where not all calls in consensus (= rounds of RT) agree?"));
  addOption(parser, seqan::ArgParseOption("m", "max-alt-qual", "Maximum value for alternative calls quality score (default value = 10)", seqan::ArgParseArgument::INTEGER, "INT"));
  
  addOption(parser, seqan::ArgParseOption("l", "max-indel-size", "Maximum allowed size for indels (default = 3, set to zero to skip all reads with indels)", seqan::ArgParseArgument::INTEGER, "INT"));

  addOption(parser, seqan::ArgParseOption("t", "threads", "How many threads to use? (default = 1)", seqan::ArgParseArgument::INTEGER, "INT"));

  addOption(parser, seqan::ArgParseOption("i", "illumina-min-score", "The character corresponding to an illumina score of zero (default: !)", seqan::ArgParseArgument::STRING, "TEXT"));
  addOption(parser, seqan::ArgParseOption("n", "no-reverse-cpt-R2", "Should the sequence of the R2 read NOT be reverse complemented?"));
  addOption(parser, seqan::ArgParseOption("s", "strand-specific", "Is the data strand specific? (observations from the two strand will not be combined if this option is set)"));
  addOption(parser, seqan::ArgParseOption("z", "zero-based", "If set, all genomic positions are reported in a zero-based numbering (one-based otherwise)"));

  addOption(parser, seqan::ArgParseOption("x", "max-mismatch-inside-consensus", "Maximum fraction of mismatched calls inside the consensus (default=0.05)", seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
  

  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  seqan::CharString bamFileName;
  seqan::getOptionValue(bamFileName, parser, "bam");

  seqan::CharString posFileName;
  seqan::getOptionValue(posFileName, parser, "pos");

  seqan::CharString outFileCalls;
  seqan::getOptionValue(outFileCalls, parser, "calls");
  seqan::CharString outFileObs;
  seqan::getOptionValue(outFileObs, parser, "observations");

  seqan::CharString outFilePrct = "consensus-quality.tab.gz";
  
  seqan::CharString r1FileName, r2FileName;
  seqan::getOptionValue(r1FileName, parser, "R1");
  if( isSet(parser, "R2") == true ){
    seqan::getOptionValue(r2FileName, parser, "R2");
  } else {
    r2FileName = "/dev/null";
  }

  seqan::CharString refGenomeFileName;
  seqan::getOptionValue(refGenomeFileName, parser, "reference-genome");

  seqan::CharString gffFileName;
  seqan::getOptionValue(gffFileName, parser, "GFF");

  int maxIndelSize = 3;
  if ( isSet(parser, "max-indel-size") == true ){
    seqan::getOptionValue(maxIndelSize, parser, "max-indel-size");
  }
  
  int minDepth = 3;
  if ( isSet(parser, "min-depth") == true ){
    seqan::getOptionValue(minDepth, parser, "min-depth");
  }
  //std::cerr<<"minDepth:"<<minDepth<<'\n';
  int minConsensusQual = 100;
  if( isSet(parser, "min-consensus-qual") == true ){
    seqan::getOptionValue(minConsensusQual, parser, "min-consensus-qual");
  }

  int nbExcludeStart = 0;
  if( isSet(parser, "exclude-start") == true ){
    seqan::getOptionValue(nbExcludeStart, parser, "exclude-start");
  }
  int nbExcludeEnd = 0;
  if( isSet(parser, "exclude-end") == true ){
    seqan::getOptionValue(nbExcludeEnd, parser, "exclude-end");
  }

  bool reverseCptSeq2 = !isSet(parser, "no-reverse-cpt-R2");

  bool allowAlternativeCalls = isSet(parser, "allow-alternative-calls");

  int maxAltQual = 10;
  if( allowAlternativeCalls == false && isSet(parser, "max-alt-qual")==true ){
    seqan::getOptionValue(maxAltQual, parser, "max-alt-qual");
  }

  char illuminaMinScore = '!';
  if( isSet(parser, "illumina-min-score")== true ){
    seqan::CharString stmp;
    seqan::getOptionValue(stmp, parser, "illumina-min-score");
    illuminaMinScore = stmp[0];
  }

  int nbThreads = 1;
  if( isSet(parser, "threads") == true ){
    seqan::getOptionValue(nbThreads, parser, "threads");
  }

  double max_csPrctMM = 0.05; // Maximum fraction of mismatches inside the consensus sequence.
  if( isSet(parser, "max-mismatch-inside-consensus") == true ){
    seqan::getOptionValue(max_csPrctMM, parser, "max-mismatch-inside-consensus");
  }
  //std::cerr<< "Maximum mismatch between calls inside consensus: "
  //<< max_csPrctMM
  //<<"\n";
      
  int oneBased = 1;
  if( isSet(parser, "zero-based") == true ){ oneBased = 0; }

  bool verbose = true;
  bool strandSpecific = false;
  int nbStrands = 1;
  if( isSet(parser, "strand-specific"  ) == true ){
    strandSpecific = true;
    nbStrands = 2;
  }
  bool combineStrands = false;
  
  int nbMaxMMperCs = 1;
  //int maxIndelSize = 6;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // END OF ARGUMENTS PARSING
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  

  // First, let's open the reference genome index:
  std::cerr<<"Opening index file from genome at '"<<refGenomeFileName<<"' ...";
  std::cerr.flush();	     
  seqan::FaiIndex faiIndex;
  if (!seqan::open(faiIndex, seqan::toCString(refGenomeFileName))){
    std::cerr << "\nERROR: Could not load FAI index for " << refGenomeFileName << "\n";
    return -1;
  }
  std::cerr<<" done.\n";
  std::cerr.flush();	     

  // Load all the positions (will need to modify this in the future to avoid loading everything in memory)
  std::cerr<<"Loading in memory all consensus positions from '"<<posFileName<<"' ...";
  std::cerr.flush();	     
  std::map<std::string, std::pair<int,int> > mPos;
  loadAllPositions(seqan::toCString(posFileName), mPos);
  std::cerr<<" done.\n";
  std::cerr.flush();	       

  // Preparing a map of all transcripts
  std::map<std::string , GFFtranscriptCoordinates> mTranscripts;
  std::cerr<<"Loading a map of all transcripts from '"<<gffFileName<<"' ...";
  std::cerr.flush();
  bool ret = loadMapOfTranscripts(seqan::toCString(gffFileName), mTranscripts);
  if( ret == false ){
    std::cerr<<"\nFAILED to load map of transcripts from '"<<gffFileName<<"'\n";
    return 1;
  }
  /*
  std::cerr<<"\n";
  mTranscripts["PCAU.43c3d.1.T00010001"].printStructure(std::cerr);
  std::cerr<<"\n";
  std::vector<int> vv;
  mTranscripts["PCAU.43c3d.1.T00010001"].getAll_tToG(vv);
  std::cerr<<"vv.size() = "<<vv.size()<<"\n";
  for(int ii=0 ; ii<vv.size() ; ii++){
    std::cerr<<"vv["<<ii<<"] -> "<<vv[ii]<<"\n";
  }
  */
  //return -2;
  std::cerr<<" done.\n";
  std::cerr.flush();



  BasesIndices BI; // This is the object that allows me to convert nucleotide into indices (and reciprocally) + reverse complement ...

  // The two objects used to store the list of calls & observations:
  CS_ScaffoldObservations cobCalls;
  CS_ScaffoldObservations cobObs;
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Preparing the gzip-compressed streams for calls and observations:
  std::ofstream ofsCalls(toCString(outFileCalls), std::ios_base::out | std::ios_base::binary); // File with all the calls made
  std::ofstream ofsObs(toCString(outFileObs), std::ios_base::out | std::ios_base::binary);  // File with all the observations made (=calls, no matter number of repeat & quality)
  std::ofstream ofsPrct(toCString(outFilePrct), std::ios_base::out | std::ios_base::binary); // File with the fraction of discordant base calls within each consensus 
  boost::iostreams::filtering_streambuf<boost::iostreams::output> outbufCalls;
  boost::iostreams::filtering_streambuf<boost::iostreams::output> outbufObs;
  boost::iostreams::filtering_streambuf<boost::iostreams::output> outbufPrct;
  outbufCalls.push(boost::iostreams::gzip_compressor());
  outbufObs.push(boost::iostreams::gzip_compressor());
  outbufPrct.push(boost::iostreams::gzip_compressor());
  outbufCalls.push(ofsCalls);
  outbufObs.push(ofsObs);
  outbufPrct.push(ofsPrct);
  std::ostream outCalls(&outbufCalls); // Stream to output the list of calls per genomic position
  std::ostream outObs(&outbufObs); // Stream to output the list of observations per genomic position
  std::ostream outPrct(&outbufPrct); // Stream to output the fraction of mismatched call inside each consensus
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  

  int BATCH_SIZE = 10000; // Not used currently (for future multi-threading version)

  /////////////////////////////////////////////////////////////////////////
  // Structures that will store the R1 & R2 read sequences + qualities:
  seqan::StringSet<seqan::CharString> ids1;
  seqan::StringSet<seqan::Dna5String> seqs1;
  seqan::StringSet<seqan::CharString> quals1;

  seqan::StringSet<seqan::CharString> ids2;
  seqan::StringSet<seqan::Dna5String> seqs2;
  seqan::StringSet<seqan::CharString> quals2;

  seqan::SeqFileIn r1FileIn(seqan::toCString(r1FileName));
  seqan::SeqFileIn r2FileIn(seqan::toCString(r2FileName));
  /////////////////////////////////////////////////////////////////////////
  

  // Now dealing with the input bam file:
  seqan::BamFileIn bamFileIn(seqan::toCString(bamFileName));
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Need to read the header before the records...
  seqan::BamHeader header;
  seqan::readHeader(header, bamFileIn);
  // Matching the numerical sequence ID to its string version:
  std::map<int32_t , std::string> mBamHeader;
  typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;
  TBamContext const & bamContext = context(bamFileIn);
  if( verbose == true ){ std::cerr<<"Processing the header from '"<<bamFileName<<"'\n"; }
  for (int i = 0; i < seqan::length(contigNames(bamContext)); ++i){
    std::string targetID(seqan::toCString(contigNames(bamContext)[i]));
    mBamHeader[i] = targetID;
  }  


  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now preparing to read the records in the BAM file
  seqan::BamAlignmentRecord record;
  //std::vector<seqan::BamAlignmentRecord> vRecords(BATCH_SIZE); // Not used currently, but maybe in a future version with multi-threading?
  int nbReadInBatch = 0;
  
  seqan::CharString ChrSeq;
  std::string previousChrID = "";
  std::string previousTranscriptID = "";
  GFFtranscriptCoordinates gt;

  seqan::CharString id1, id2;
  seqan::CharString seq1, seq2;
  seqan::CharString qual1, qual2;

  std::vector<int> vTtoG; // Vector that will store the transcript to genomic coordinates mappin information.

  // Reading bam records one at a time...
  while (!seqan::atEnd(bamFileIn)){
    seqan::readRecord(record, bamFileIn);

    bool mappedRC = seqan::hasFlagRC(record);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Extracting information from the TAGs
    
    // Extracting the chromosome name from the tag ("MS")
    seqan::BamTagsDict tagsDict(record.tags);
    unsigned tagIdx = 0;
    if( !seqan::findTagKey(tagIdx, tagsDict, "MS") ){
        std::cerr << "ERROR: entry does not contain chromosome information!\n";
	return -2;
    }
    seqan::CharString chromTag;
    if( !seqan::extractTagValue(chromTag, tagsDict, tagIdx) ){
        std::cerr << "ERROR: There was an error extracting chromosome from tags!\n";
	return -2;
    }
    std::string chrID(seqan::toCString(chromTag));

    // Extracting the breakpoint location ("BP"):
    if( !seqan::findTagKey(tagIdx, tagsDict, "BP") ){
      std::cerr << "ERROR: entry does not contain breakpoint information!\n";
      return -2;
    }
    int posBreakPoint;
    if( !seqan::extractTagValue(posBreakPoint, tagsDict, tagIdx) ){
        std::cerr << "ERROR: There was an error extracting the breakpoint position from tags!\n";
	return -2;
    }

    // Transcript's strand (ST):
    if( !seqan::findTagKey(tagIdx, tagsDict, "ST") ){
      std::cerr << "ERROR: entry does not contain transcript's strand information!\n";
      return -2;
    }
    char tStrand;
    if( !seqan::extractTagValue(tStrand, tagsDict, tagIdx) ){
        std::cerr << "ERROR: There was an error extracting the transcript's strand information from tags!\n";
	return -2;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // If the its a new chromosome, I need to print the previous one and read the sequence of the new one + prepare a new CS_ScaffoldObservations object.
    if( chrID != previousChrID ){
      std::cerr<<"New chromosome: '"<<chrID<<"'\n";
      if( previousChrID != "" ){
	cobCalls.print(outCalls, true, false, false, combineStrands, ChrSeq, BI, oneBased);
	cobObs.print(outObs, true, false, false, combineStrands, ChrSeq, BI, oneBased);
      }
      previousChrID = chrID;
      unsigned idx = 0;
      if( !getIdByName(idx, faiIndex, seqan::toCString(chrID.c_str())) ){
	std::cerr << "ERROR: FAI index has no entry for '"<<chrID<<"' --> ABORTING PROGRAM.\n";
	return -1;
      }
      int seqLength = sequenceLength(faiIndex, idx);
      cobCalls.init(seqLength, chrID, nbStrands);
      cobObs.init(seqLength, chrID, nbStrands);
      seqan::clear(ChrSeq);
      readSequence(ChrSeq, faiIndex, idx);
    } // Done preparing for the new chromosome...
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // The REF field of the bam entry should be a transcript ID (consensus sequences were mapped against transcripts)
    std::string transcriptID = mBamHeader[record.rID];
    //std::string sTmp;
    if( transcriptID != previousTranscriptID ){
      // New transcript
      previousTranscriptID = transcriptID;
      gt = mTranscripts[transcriptID];
      vTtoG.clear();
      bool res_tToG = gt.getAll_tToG(vTtoG);

      //std::stringstream issTmp;
      //issTmp << record.qName;
      //issTmp >> sTmp;
      
      if( res_tToG == false ){
	std::cerr<<"PROBLEM DURING CONVERSTION OF TRANSCRIPT to GENOMIC COORDINATES WITH: "<<record.qName<<"\n";
	std::cerr<<"Printing the structure for "<<transcriptID<<"\n";
	mTranscripts[transcriptID].printStructure(std::cerr);
	std::flush(std::cerr);
	return -1;
	//continue;
      }
    }


    // Read the corresponding read sequence (the fastq files with read/pair must be exactly in the same order as the bam file)
    std::string readID(toCString(record.qName));
    
    //std::cerr<<readID<<'\n';
    //std::cerr.flush();
    readRecord(id1, seq1, qual1, r1FileIn);
    if( pairedFile == true ){
      readRecord(id2, seq2, qual2, r2FileIn);
    }
    if( pairedFile==true && (id1 != id2 || id1 != record.qName) ){
      std::cerr<<"ERROR: BAM file and reads file are out of sync.\nID in bam: '"<<record.qName<<"'\nID in fastq: '"<<id1<<"' and '"<<id2<<"'\n";
      return 1;
    }

    if( seqan::length(record.cigar) > 1 ){

      // I work only on reads with a cigar string that has th form: [numeric]M[numeric]D/I[numeric]M
      if( seqan::length(record.cigar) != 3 ||
	  record.cigar[0].operation != 'M' ||
	  record.cigar[2].operation != 'M' ||
	  (record.cigar[1].operation != 'I' && record.cigar[1].operation != 'D') ||
	  record.cigar[1].count > maxIndelSize
	  )
	continue;

      /**/
      if( seqan::length(record.cigar) == 3 && record.cigar[0].operation=='M' && record.cigar[2].operation=='M' ){
	// Extracting the value in the middle of the cigar string
	char rType = record.cigar[1].operation;
	int indelSize = record.cigar[1].count;
	int indelPosInCs = record.cigar[0].count;
	int indelCsLg = seqan::length(record.seq);
	if( rType == 'D' || rType == 'I' ){
	  int indelGenomicPos = vTtoG[record.beginPos + indelPosInCs]; // Check if this this should be rPos+1 or maybe rPos-1
	  string indelType = "INS";
	  if( rType == 'D' ){ indelType = "DEL"; }
	  cout<<indelType<<'\t'
	      <<record.qName<<'\t'
	      <<indelPosInCs<<'\t'
	      <<indelCsLg<<'\t'
	      <<chrID<<'\t'
	      <<indelGenomicPos<<'\t'
	      <<"B_FROM"<<'\t'
	      <<indelSize<<'\t'
	      <<tStrand<<'\t'
	      <<mappedRC
	      <<'\n';
	}
      }      
      //std::cerr<<"SKIPPING "<<record.qName<<'\n';
      // WARNING: WITH THIS VERSION, I DO NOT CALL BASE-SUBS IN READS WITH INDELS!
      continue;
      /**/
    }    
    
    // Now it is time to actually process the read
    //std::cerr<<"Processing: "<<record.qName<<"\n";

    int csLg = seqan::length(record.seq);
    std::pair<int,int> ppos = mPos[readID];
    std::vector<std::string> vSeq(2);
    std::vector<std::string> vQual(2);
    std::vector<int> vPosFirstBreakPoint(2);

    int posBreakPoint1 = posBreakPoint;
    if( mappedRC == true ){
      seqan::reverseComplement(seq1);
      seqan::reverse(qual1);
      posBreakPoint1 = ( length(seq1) % csLg ) + posBreakPoint;
    }

    int posBreakPoint2 = (ppos.second + posBreakPoint) % csLg;
    if( pairedFile == true ){
      if( (reverseCptSeq2 == true && mappedRC == false) || (reverseCptSeq2==false && mappedRC==true)  ){
	//std::cerr<<"REVERSE!\n";
	seqan::reverseComplement(seq2);
	seqan::reverse(qual2);
      }

      if( mappedRC == true ){
	posBreakPoint2 = ( (length(seq2) - ppos.second) + posBreakPoint) % csLg;
      }
      
      vSeq[1] = std::string(seqan::toCString(seq2));
      vQual[1] = std::string(seqan::toCString(qual2));
      vPosFirstBreakPoint[1] = posBreakPoint2;
    }
      
    //std::cerr<<"seq2:'"<<seq2<<"'\n";
    vSeq[0] = std::string(seqan::toCString(seq1));
    vQual[0] = std::string(seqan::toCString(qual1));
    vPosFirstBreakPoint[0] = posBreakPoint1;

    /*
    std::cerr<<"seq1:"<<seq1<<'\n';
    std::cerr<<"seq2:"<<seq2<<'\n';
    std::cerr<<"ppos: "    <<ppos.first<<" - "<<ppos.second<<'\n';
    std::cerr<<"BP1:"<<posBreakPoint<<'\n';
    std::cerr<<"BP2:"<<posBreakPoint2<<'\n';
    */
    //auto start1 = high_resolution_clock::now();
    CS_Consensus cs(vSeq, vQual, csLg, vPosFirstBreakPoint, BI, illuminaMinScore);
    int nbMatch, nbMismatch;
    double csPrctMM = cs.getPrctMismatchesInsideConsensus(&nbMatch, &nbMismatch);
    outPrct<<readID<<'\t'<<csPrctMM<<"\n";
    if( csPrctMM > max_csPrctMM )
      continue;

    //auto stop1 = high_resolution_clock::now();
    //auto duration1 = duration_cast<nanoseconds>(stop1 - start1);
    //std::cerr<<"Building the consensus: "<<duration1.count()<<'\n';

    //auto start2 = high_resolution_clock::now();
    int nbPosDiff = callThisConsensus(cs,
				      nbExcludeStart, nbExcludeEnd,
				      record,
				      chrID, tStrand, mappedRC, strandSpecific,
				      ChrSeq, gt, vTtoG,
				      cobCalls, cobObs,
				      minDepth, allowAlternativeCalls, maxAltQual, minConsensusQual, nbMaxMMperCs,
				      BI, oneBased);

    if( nbPosDiff < 0 ){
      std::cerr<<"ERROR: nbPosDiff<0 --> Printing the structure for '"<<transcriptID<<"'\n";
      mTranscripts[transcriptID].printStructure(std::cerr);
      std::flush(std::cerr);
      exit(1);
    }
    
    //auto stop2 = high_resolution_clock::now();
    //auto duration2 = duration_cast<nanoseconds>(stop2 - start2);
    //std::cerr<<"Making the calls: "<<duration2.count()<<'\n';


    //std::cerr<<readID<<"\t"<<nbPosDiff<<'\n';
    //cs.printDetails(std::cout, BI);
    //std::cerr<<cs.getConsensus()<<'\n';   
  }
  
  cobCalls.print(outCalls, true, false, false, combineStrands, ChrSeq, BI, oneBased);
  cobObs.print(outObs, true, false, false, combineStrands, ChrSeq, BI, oneBased);

  boost::iostreams::close(outbufCalls);
  boost::iostreams::close(outbufObs);
  boost::iostreams::close(outbufPrct);
  ofsCalls.close();
  ofsObs.close();
  ofsPrct.close();
  
  return 0;
}



int callThisConsensus(CS_Consensus &cs,
		      int nbExcludeStart, int nbExcludeEnd,
		      seqan::BamAlignmentRecord &record,
		      std::string &chrID, char tStrand, bool mappedRC, bool strandSpecific,
		      seqan::CharString &ChrSeq,
		      GFFtranscriptCoordinates &gt,
		      std::vector<int> &vTtoG,
		      CS_ScaffoldObservations &cobCalls, CS_ScaffoldObservations &cobObs,
		      int minDepth, int allowAlternativeCalls, int maxAltQual, int minConsensusQual,
		      int nbMaxMMperCs, // How many divergent base calls are allowed per consensus?
		      BasesIndices &BI,
		      int oneBased
		      ){


  std::string readID(toCString(record.qName));
  //std::cerr<<"minDepth in call:"<<minDepth<<'\n';  
  int nbPosCallable = 0;
  int nbPosDiff = 0;
  int nbLocalCandidates = 0;
  int csLg = cs.getLength();
  
  std::vector<int> vPosCallable(csLg);
  //std::vector<int> vPosCandidates(csLg);
  std::map<int, CS_Candidate> mCand;
  
  std::map<int, char> mCalls;
  std::map<int, char> mObs;

  char gMappedStrand = '+';
  bool reverseGenome = false;
  bool reverseCs = false;

  if( strandSpecific == false ){
    if( tStrand == '-' ){
      reverseCs = true;      
    }
  } else {
    if( tStrand == '+' ){
      if( mappedRC == false ){
	// consensus maps on forward strand
	;	
      } else {
	// consensus maps on reverse strand
	reverseGenome = true;
	reverseCs = true;
	gMappedStrand = '-';
      }
    } else {
      // For transcripts on the reverse strand:
      if( mappedRC == false ){
	reverseGenome = true;
	gMappedStrand = '-';
      } else {
	// read actually maps on the forward strand
	reverseCs = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Doing a first round to check that the consensus does not have too many mismatches
  // (and storing the results in temporary variables to flush into the observations/candidates if the consensus passes the threshold for max number of mismatches

  int cPos = 0;
  int indelOffset = 0;
  std::vector<bool> vStatus_csPositions(csLg, false); // This vector will record the list of positions (in the consensus) that pass all the thresholds (initialized at false for every position)
  
  for(int iic = 0 ; iic < seqan::length(record.cigar) ; iic++){
    if( record.cigar[iic].operation == 'M' ){
      int matchLength = record.cigar[iic].count;
  
      for(int ii = 0 ; ii < matchLength ; ii++, cPos++){

	if( cPos >= csLg ){
	  std::cerr<<"ERROR: cPos ("
		   <<cPos
		   <<") >= csLg ("
		   <<csLg
		   <<") for: "<< record.qName
		   <<"\n";
	  exit(1);
	}
	
	int tPos = record.beginPos + cPos + indelOffset;
	int gPos = vTtoG[tPos];

	if( gPos > seqan::length(ChrSeq) ){
	  std::cerr<<"ERROR OF LENGTH WITH "<<record.qName<<"\n";
	  std::cerr<<"vTtoG["<<tPos<<"] = "<<gPos<<"\n";
	  std::flush(cerr);
	  return -1;
	}
	char gBase = ChrSeq[gPos];
	if( reverseGenome == true ){ gBase = BI.getBaseCpt(gBase); }
	//std::cerr<<cPos<<" -> "<<tPos<<" -> "<<gPos<<" -> "<<record.seq[cPos]<<" - "<<gBase<<'\n';
	int depth = -1;
	char csBase = 'N';
	int csQual = -1;
	int altQual = -1;
	bool resGetParam = cs.getParameters(cPos, &csBase, &csQual, &altQual, &depth);
	if( resGetParam == false ){
	  std::cerr<<"ERROR WHILE GETTING PARAMTERS OF POSITION IN CONSENSUS.\n";
	  return -1;
	}
	if( reverseCs == true ){ csBase = BI.getBaseCpt(csBase); }
	mObs[gPos] = csBase;

    
	if( csBase != gBase ){ nbPosDiff++; }
	if( cPos>nbExcludeStart && cPos<(csLg-nbExcludeEnd) && depth>=minDepth && (altQual==0 || allowAlternativeCalls==true) && altQual<=maxAltQual && (csQual-altQual)>=minConsensusQual ){
	  vPosCallable[nbPosCallable++] = cPos;
	  vStatus_csPositions[cPos] = true;
	  mCalls[gPos] = csBase;
	  if( csBase != gBase ){
	    CS_Candidate cd(gPos, cPos, gBase, csBase);
	    mCand[cPos] = cd;
	    nbLocalCandidates++;
	  }
	}
      }
    }
    // If this is an insertion, I need to move the position on the consensus (=cPos) of the size of the insertion and reduce the offsetIndel by the same amount.
    // For example, if the cigar string is: 23M2I61M I should reach this point when cPos = 23. The next two nucleotides in the consensus will not match anything in the genome (insertion), so the next
    // cPos I need to look t is (23+2)=25. However, the position in the transcript is still 23 (--> this is why I use the indelOffset variable).
    if( record.cigar[iic].operation == 'I' ){
      indelOffset -= record.cigar[iic].count;
      cPos += record.cigar[iic].count;
    }
    if( record.cigar[iic].operation == 'D' ){
      indelOffset += record.cigar[iic].count;
    }
    
  }

  if( seqan::length(record.cigar) == 3 ){
    // If there is an indel in this consensus, I will check that the nucleotides around it fulfill all the criteria (depth, ...)
    // TO DO !!!
    ;
  }
  
  if( nbPosDiff > MAX_POS_DIFF_BEFORE_WARNING ){
    std::cerr<<"WARNING: '"<<readID<<"' has "<<nbPosDiff<<" mismatches with the reference genome!\n";
  }

  int iStrand = BI.getStrandIndice(gMappedStrand);
  std::map<int,char>::iterator itm;

  // Adding the validated observations (= calls)
  if( nbPosDiff <= nbMaxMMperCs ){
    for(itm=mCalls.begin() ; itm!=mCalls.end() ; itm++){
      //std::cerr<<"addCall("<<itm->second<<","<<itm->first<<","<<iStrand<<")\n";
      cobCalls.addObservation(itm->second, 1, itm->first, iStrand, BI);
    }
  }
  // Adding all the observations (including low quality)
  for(itm=mObs.begin() ; itm!=mObs.end() ; itm++){
    cobObs.addObservation(itm->second, 1, itm->first, iStrand, BI);
  }

  // Printing out the candidates!
  std::string errorType("BASE_SUB");
  std::ostream &outCand = std::cout;
  if( nbPosDiff <= nbMaxMMperCs && nbLocalCandidates > 0 ){
    std::map<int, CS_Candidate>::iterator itc;
    for(itc=mCand.begin() ; itc!=mCand.end() ; itc++){
      (itc->second).print(readID, chrID, csLg, outCand, tStrand, mappedRC, oneBased);
    }
  }

  return nbPosDiff;
}
