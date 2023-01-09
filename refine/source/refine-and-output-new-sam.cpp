//###############################################################################################
//
// This is the program that parses the bam output from Kallisto and finds the correct breakpoint


#include "unistd.h"

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/parallel.h>

#include <edlib.h>

#include <map>
#include <vector>
#include <iterator>

#include <omp.h>

#include <gff-utils.hpp>

#define NB_READS_PER_PACKET 10000 // 10000
#define NB_MAX_ALN_PER_READ 200 // 2000

using namespace seqan;



// rc_findMinEdit
//! Finds the best alignment between a read and its target.
/*!
This function will try all possible breakpoint position in a given read and return the best alignment.
\param vBam a vector of bamAlignmentRecord (from seqan) that contains all the bam alignments for a given read
\param nbAln The number of alignments to consider in the vector (because the vectors are pre-allocated for faster memory usage, I keep track of the number of elements in it "by hand")
\param mRefSeqs A map (STL) containing the DNA sequence of each possible target (the key is an integer corresponding to the number in the bam header)
\param posOfBestAln Will be set to the position of the breakpoint corresponding to the best alignment.
 */
int rc_findMinEdit(
		   std::vector<seqan::BamAlignmentRecord> &vBam,
		   int nbAln,
		   std::map<int32_t , std::pair<std::string , std::string> > &mRefSeqs,
		   //std::map<int32_t , const char *> &mRefSeqs,
		   int *posOfBestAln,
		   int *startInTranscript,
		   int *whichRecordOfBestAln,
		   std::string &transcriptID,
		   std::string &debreakedReadSeq,
		   std::string &debreakedReadQual,
		   std::string &cigar,
		   bool printDetails
		   );



bool fillCigarElements(seqan::BamAlignmentRecord &recordToOutput, std::string &cigar);

int main(int argc, char const *argv[]){

  int i, j;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Argument parsing.
  seqan::ArgumentParser parser("refine");

  //addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption(
					  "k", "kallisto-bam", "The sam/bam file produced by Kallisto",
					  seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption(
					  "r", "reference", "The fasta reference file used for mapping",
					  seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption(
					  "g", "gff", "The GFF/GTF file containing the location of all the transcripts in the reference file.",
					  seqan::ArgParseArgument::STRING, "TEXT"));

  addOption(parser, seqan::ArgParseOption(
					  "t", "threads", "Number of threads to use.",
					  seqan::ArgParseArgument::INTEGER, "INT"));


  addOption(parser, seqan::ArgParseOption("p", "print-details", "For debugging use only."));

  
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;
  
  seqan::CharString bamFileName;
  seqan::getOptionValue(bamFileName, parser, "kallisto-bam");

  seqan::CharString refFile;
  seqan::getOptionValue(refFile, parser, "reference");

  seqan::CharString gffFile;
  seqan::getOptionValue(gffFile, parser, "gff");

  int NB_THREADS;
  seqan::getOptionValue(NB_THREADS, parser, "threads");

  bool printDetails = false;
  if( isSet(parser, "print-details") == true ){ printDetails = true; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Storing the correspondance between transcriptID and chromosome/scaffold from the GFF/GTF file
  // This will be used later on to add the scaffold/chromosome information to the new sam/bam file.
  // Doing so will then allow to sort the bam file by chromosome
  std::map<std::string, std::pair<std::string , char> > mTlink;
  fprintf(stderr, "Pairing each transcript to its chromosome and strand ...");
  getScaffoldAndStrandPerTranscript(mTlink, seqan::toCString(gffFile));
  fprintf(stderr, " done.\n");
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Reading the fasta reference (this should be a fasta file containing the sequence of every transcript):
  seqan::StringSet<seqan::CharString> refIDs;
  seqan::StringSet<seqan::Dna5String> refSeqs;
  fprintf(stderr, "Reading sequence from the reference fasta: %s ...", seqan::toCString(refFile));
  seqan::SeqFileIn seqFileIn(seqan::toCString(refFile));
  seqan::readRecords(refIDs, refSeqs, seqFileIn);
  fprintf(stderr, " done.\n");


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
  fprintf(stderr, "Processing the header from '%s' ...", seqan::toCString(bamFileName));
  std::map<std::string , int32_t > mHeader;
  // Saving the correspondance between sequence number in the bam header and corresponding ID to access in the fasta file:
  for (i = 0; i < seqan::length(contigNames(bamContext)); ++i){
    std::string targetID(seqan::toCString(contigNames(bamContext)[i]));
    mHeader[targetID] = i;
    //std::cout << contigNames(bamContext)[i] << '\t'
    //    <<"'"<<targetID<<"'" << '\t'
    //    << contigLengths(bamContext)[i] << '\n';
  }
  fprintf(stderr, " done.\n");


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Linking refSeqs to their corresponding bam header int ID
  fprintf(stderr, "Loadding reference sequences in map with bam header correspondance ...");
  std::map< int32_t , std::pair<std::string , std::string> > mRefSeqs;
  //std::map< int32_t , const char * > mRefSeqs;
  for(i = 0 ; i < length(refIDs) ; i++){
    std::string refSeqID(seqan::toCString(refIDs[i]));
    std::stringstream ss;
    ss << refSeqs[i];
    std::string refSeqString = ss.str();
    std::map< std::string , int32_t >::iterator itmh = mHeader.find(refSeqID);
    if( itmh == mHeader.end() ){
      fprintf(stderr, "\nERROR: '%s' not found in the bam header --> EXITING!\n", refSeqID.c_str());
      return -2;
    }
    int32_t seqID_int = itmh->second;
    std::pair<std::string, std::string> ppr(refSeqID, refSeqString);
    mRefSeqs[seqID_int] = ppr;
  }
  fprintf(stderr, " done.\n");

  fprintf(stderr, "Clearing the stringsets ...");
  seqan::clear(refIDs);
  seqan::clear(refSeqs);
  fprintf(stderr, " done.\n");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Preparing the (bam) output file:
  seqan::BamFileOut bamFileOut(seqan::context(bamFileIn), std::cout, seqan::Bam());
  seqan::writeHeader(bamFileOut, header);
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Read records.
  int nbTotalBamRecordsRead = 0;
  int pos;
  bool atEndOfBamFile = false;
  seqan::BamAlignmentRecord record;

  fprintf(stderr, "Allocating memory to store a packet of BamAlignmentRecords ...");
  std::vector< std::pair<int, std::vector<seqan::BamAlignmentRecord> > > vRec(NB_READS_PER_PACKET+1);
  std::vector< std::pair<int, std::vector<seqan::BamAlignmentRecord> > >::iterator itr;
  std::vector< seqan::BamAlignmentRecord > vBamOutput(NB_READS_PER_PACKET+1);
  // Initializing the vector that will store the packets of records (records are kept in memory as large packets to be then processed in parallel)
  for(i=0 ; i<= NB_READS_PER_PACKET ; i++){
    std::vector< seqan::BamAlignmentRecord > vtmp(NB_MAX_ALN_PER_READ);
    std::pair< int , std::vector< seqan::BamAlignmentRecord > > pp(0, vtmp);
    vRec[i] = pp;
  }
  fprintf(stderr, " done.\n");
  
  std::string lastWarningReadID = "";
  //std::cerr<<"REFINE\n";
  // Start of the main loop that reads the bam file one entry at a time.
  while( atEndOfBamFile == false ){
    int nbBamReadsRead = 0;
    std::string previousID = "";

    // Start of the inner loop that reads one packet of bam records.
    while( atEndOfBamFile == false && nbBamReadsRead < NB_READS_PER_PACKET ){
      
      readRecord(record, bamFileIn);
      //std::cerr<<"readID:'"<<record.qName<<"'\n";

      nbTotalBamRecordsRead++;
      if( nbTotalBamRecordsRead % 1000000 == 0 ){
	fprintf(stderr, "%d bam records read ...\n", nbTotalBamRecordsRead);
      }
      
      std::string readID(seqan::toCString(record.qName));
      if( readID != previousID ){
	itr = vRec.begin() + nbBamReadsRead;
	nbBamReadsRead++;
	//std::cerr<<"readID:'"<<readID<<"' (nbBamReadsRead:"<<nbBamReadsRead<<")\n";
	previousID = readID;
	
	if( nbBamReadsRead == NB_READS_PER_PACKET ){
	  itr = vRec.end();
	  //std::cerr<<"WARNING, LAST READ IN THE PACKET!!!\n";
	  // DO SOMETHING --> SAVE bamAlignmentRecord into a temporary variable which will be incorporated at position 1 of the structure for the next iteration of the big loop
	} else {
	  (*itr).first = 0;
	}
      }

      if( itr != vRec.end() ){
	// Add the bamAlignmentRecord to the structure:
	pos = (*itr).first;
	if( pos < NB_MAX_ALN_PER_READ ){
	  ((*itr).second)[pos] = record;
	  (*itr).first = pos + 1;
	} else {
	  if( readID != lastWarningReadID ){
	    fprintf(stderr, "WARNING: MORE THAN %d ALN FOR %s\n", NB_MAX_ALN_PER_READ, readID.c_str());
	    lastWarningReadID = readID;
	  }
	}
	
      }

      atEndOfBamFile = seqan::atEnd(bamFileIn);

    }

    int nbBamRecordsToOutput = 0;
    // Process the packet of bam records
    //fprintf(stderr, "Processing the packet ...\n");
    bool pbStop = false;
    #pragma omp parallel num_threads(NB_THREADS)
    #pragma omp for
    for(int iRecord = 0 ; iRecord<nbBamReadsRead ; iRecord++){
      //#pragma omp cancellation point parallel
      int posOfMinEdit = -1;
      int startInTranscript = -1;
      int whichRecordOfMinEdit = -1;
      std::string debreakedReadSeq;
      std::string debreakedReadQual;
      std::string transcriptID;
      std::string cigar = "";


      // This is the function that will try all possible arrangements of the circular repeat to find the best possible alignment to the transcript.
      int minEdit = rc_findMinEdit(vRec[iRecord].second, vRec[iRecord].first, mRefSeqs, &posOfMinEdit, &startInTranscript, &whichRecordOfMinEdit, transcriptID, debreakedReadSeq, debreakedReadQual, cigar, printDetails);
      std::string readID( seqan::toCString(vRec[iRecord].second[0].qName) );
      //std::cerr<<"readID:'"<<readID<<" --> "<<minEdit<<"'\n";

      // If a good alignment was found, then it needs to be outputed in bam
      if( posOfMinEdit >= 0 ){
	seqan::BamAlignmentRecord recordToOutput = vRec[iRecord].second[whichRecordOfMinEdit];
	//std::cerr<<"ReadID: '"<<recordToOutput.qName<<"' (iRecord="<<iRecord<<") - whichRecordofMinEdit = "<<whichRecordOfMinEdit<<"\n"; // <- readID is now empty for rRNAs!!!
	seqan::clear(recordToOutput.seq);
	seqan::clear(recordToOutput.qual);
	recordToOutput.seq = debreakedReadSeq;
	recordToOutput.qual = debreakedReadQual;
	//recordToOutput.tags += "MS:Z:12";
	seqan::BamTagsDict tagsDict(recordToOutput.tags);
	seqan::eraseTag(tagsDict, "ZW");
	seqan::setTagValue(tagsDict, "NH", minEdit);
	seqan::setTagValue(tagsDict, "BP", posOfMinEdit);
	//seqan::setTagValue(tagsDict, "TS", startInTranscript);
	      
	std::map<std::string , std::pair<std::string, char> >::iterator itl = mTlink.find(transcriptID);
	if( itl == mTlink.end() ){
          //#pragma omp critical
	  fprintf(stderr, "ERROR: '%s' not found in map of transcripts location/strand.\n", transcriptID.c_str());
	  //std::cerr<<"Map dump:\n";
	  //for(itl=mTlink.begin() ; itl!=mTlink.end() ; itl++){
	  //fprintf(stderr, "'%s'\t->\t['%s' , '%c']\n", (itl->first).c_str(), (itl->second).first.c_str(), (itl->second).second);
	  //}
	  pbStop = true;
	  //	  #pragma omp cancel parallel
	} else {
	  std::string chrID = (itl->second).first;
	  char strand = (itl->second).second;
	  seqan::setTagValue(tagsDict, "ST", strand);	
	  seqan::setTagValue(tagsDict, "MS", chrID.c_str());
	  fillCigarElements(recordToOutput, cigar);
	  recordToOutput.beginPos = startInTranscript;
	  //std::cerr<<"StartInTranscript: "<<startInTranscript<<'\n';
	  //vBamOutput[nbBamRecordsToOutput] = recordToOutput;
	  //nbBamRecordsToOutput++;
#pragma omp critical
	  {
	    seqan::writeRecord(bamFileOut, recordToOutput);
	  }
	  //std::cerr<<"tags:\t"<<recordToOutput.tags<<"\n";
	}
      }
      //fprintf(stderr, "%d\t%s\t%d\t%d\t%s\n",vRec[i].first, readID.c_str(), minEdit, posOfMinEdit, cigar.c_str());
    }

    // MODIFY THIS TO ABORT THE PROGRAM IN CASE A TRANSCRIPT IS NOT FOUND.
    if( pbStop == true ){
      //return -1;
      ;
    }

    //for(i=0 ; i<nbBamRecordsToOutput ; i++){
    //seqan::writeRecord(bamFileOut, vBamOutput[i]);
    //seqan::clear(vBamOutput[i]);
    //}
    
    //fprintf(stderr, "Clearing the structure...");
    // Clear the structure that stores the bam records:
    for(i=0 ; i<= NB_READS_PER_PACKET ; i++){
      int nbToClear = vRec[i].first;
      for(j=0 ; j<nbToClear ; j++){
	seqan::clear(vRec[i].second[j]);
      }
      vRec[i].first = 0;
    }
    //fprintf(stderr, "done.\n");
    
  }
  

  return 0;
}


int rc_findMinEdit(std::vector<seqan::BamAlignmentRecord> &vBam,
		   int nbAln,
		   //std::map<int32_t , const char *> &mRefSeqs,
		   std::map<int32_t , std::pair<std::string , std::string> > &mRefSeqs,
		   //seqan::StringSet<seqan::CharString> & refIDs, 
		   //seqan::StringSet<seqan::CharString> & refSeqs,
		   int *posOfBestAln, // position in the consensus where the best alignment starts (= breakpoint)
		   int *startInTranscript,
		   int *whichRecordOfBestAln,
		   std::string &transcriptID,
		   std::string &debreakedReadSeq,
		   std::string &debreakedReadQual,
		   std::string &cigar,
		   bool printDetails
		   ){

  int EXTRACT_PADDING = 100;
  int recordNum, posStart;
  int bestRecord = -1;
  int minDist = -1;
  int bestBreakPoint = -1;
  int bestStartInTranscript = -1;
  int nb = vBam.size();
  int EDLIB_MAX_EDIT_DISTANCE_SEARCH = 2;
  
  bool done = false;
  
  std::stringstream ss;
  ss << vBam[0].seq;
  std::string doubleReadSeqString = ss.str();
  int readLength = doubleReadSeqString.size() / 2;

  std::string readID(seqan::toCString(vBam[0].qName));
  
  //std::map<int32_t, const char *>::iterator itm;
  std::map<int32_t, std::pair<std::string , std::string> >::iterator itm;
  for(recordNum=0 ; !done && recordNum<nbAln ; recordNum++){
    int32_t rID = vBam[recordNum].rID;
    itm = mRefSeqs.find(rID);
    if( itm == mRefSeqs.end() ){
      //fprintf(stderr, "WARNING: rID %d not found (%s).\n", rID, toCString(vBam[recordNum].qName));
      return -2; // <- This indicates that the corresponding read did not have a hit in the sam/bam file to start with.
    }

    // Extract the sub-region of the transcript where the alignment should take place.
    int mapStart = vBam[recordNum].beginPos;
    int alnLength = getAlignmentLengthInRef(vBam[recordNum]);
    int posStartExtract = (mapStart - readLength) - EXTRACT_PADDING;
    int actualLeftPadding = EXTRACT_PADDING;
    if( posStartExtract < 0 ){
      actualLeftPadding = EXTRACT_PADDING + posStartExtract;
      //std::cerr<<"Actual padding:"<<actualLeftPadding<<'\n';
      posStartExtract = 0;
    }
    int posEndExtract = mapStart + alnLength + readLength + EXTRACT_PADDING;
    int refSeqLength = (itm->second).second.size();
    if( posEndExtract >= refSeqLength ){ posEndExtract = refSeqLength-1; }
    int extractLength = (posEndExtract - posStartExtract);
    std::string tSeqString = (itm->second).second.substr(posStartExtract, extractLength);
    if( printDetails == true ){
      std::cerr<<"Aligning against: \n"<<tSeqString<<'\n';
    }
    const char *tSeq =  tSeqString.c_str(); //(itm->second).substr(posStart, extractLength).c_str();   
    int tLength = strlen(tSeq);
    
    // Now, perform the alignments:
    // Generating all possible rotations (= all possible break points)
    for(posStart = 0 ; !done && posStart < readLength ; posStart++){
      //fprintf(stderr, "posStart:%d - readLength:%d\n", posStart, readLength);
      std::string readSeqString = doubleReadSeqString.substr(posStart, readLength);
      const char *readSeq = readSeqString.c_str();
      
      EdlibAlignResult result = edlibAlign(readSeq, readLength, tSeq, tLength,
      	   edlibNewAlignConfig(EDLIB_MAX_EDIT_DISTANCE_SEARCH, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

      //EdlibAlignResult result = edlibAlign(readSeq, readLength, tSeq, tLength,
      //		   edlibNewAlignConfig(EDLIB_MAX_EDIT_DISTANCE_SEARCH, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

      //fprintf(stderr, "Query: %s\nTarget:%s\n", readSeq, tSeq);
      
      if (result.status == EDLIB_STATUS_OK){
	if( result.editDistance >= 0 ){
	  if( minDist<0 || result.editDistance < minDist ){
	    minDist = result.editDistance;
	    bestBreakPoint = posStart;
	    bestRecord = recordNum;
	    transcriptID.clear();
	    transcriptID = (itm->second).first;
	    
	    edlibFreeAlignResult(result);
	    EdlibAlignResult resultPath = edlibAlign(readSeq, readLength, tSeq, tLength,
	    				     edlibNewAlignConfig(EDLIB_MAX_EDIT_DISTANCE_SEARCH, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	    if( resultPath.status != EDLIB_STATUS_OK ){
	      std::cerr<<"ERROR WHILE COMPUTING THE PATH ALIGNMENT!\n";
	      return -1;
	    }
	    const char* c_cigar = edlibAlignmentToCigar(resultPath.alignment, resultPath.alignmentLength, EDLIB_CIGAR_STANDARD);
	    //std::cerr<<"result.editDistance:"<<resultPath.editDistance<<" - Number of locations: "<<resultPath.numLocations<<'\n';
	    //std::cerr<<resultPath.startLocations[0]<<'\n';
	    bestStartInTranscript = resultPath.startLocations[0] + posStartExtract;

	    if( printDetails == true ){
	      std::cerr<<"New best alignment (pos in DCS: "<<posStart
		       <<" - start aln in extracted tSeq: "<<resultPath.startLocations[0]
		       <<" - start in transcript:"<<bestStartInTranscript
		       <<" - edit distance:"<<minDist
		       <<") with: "<<readSeq<<'\n';
	    }
	    
	    cigar.clear();
	    cigar.assign(c_cigar);
	    edlibFreeAlignResult(resultPath);
	    if( minDist == 0 ){ done = true; }
	  }
	}
      }
      
      //edlibFreeAlignResult(result);

    } // Next breakpoint position
    
  } // Next BamAlignmentRecord

  if( minDist >= 0 ){
    // Passing all the results back to the program that called the function
    //(*posOfBestAln) = posStart;
    (*posOfBestAln) = bestBreakPoint; // <- Might need to change this?
    (*startInTranscript) = bestStartInTranscript;
    (*whichRecordOfBestAln) = bestRecord;
    std::string readSeq = doubleReadSeqString.substr(bestBreakPoint, readLength);
    std::stringstream ssq;
    ssq << vBam[0].qual;
    std::string doubleReadQualString = ssq.str();
    std::string readQual = doubleReadQualString.substr(bestBreakPoint, readLength);
    debreakedReadSeq = readSeq;
    debreakedReadQual = readQual;

    //char *c_cigar = "39M";
 
  }
  
  return minDist;
}



bool fillCigarElements(seqan::BamAlignmentRecord &recordToOutput, std::string &cigar){

  seqan::clear(recordToOutput.cigar);
  
  int lg = cigar.size();
  for(int i=0 ; i<lg ; i++){
    std::string sCount = "";
    char c = cigar[i];
    while( c>='0' && c<='9' ){
      sCount.push_back(c);
      i++;
      c = cigar[i];
    }

    int nb = atoi(sCount.c_str());
    seqan::CigarElement<> element;
    element.count = nb;
    element.operation = c;
    seqan::appendValue(recordToOutput.cigar, element);    
    
  }

  return true;
}
