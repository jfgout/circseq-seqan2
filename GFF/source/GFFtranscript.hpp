#ifndef _GFF_TRANSCRIPT_HPP_
#define _GFF_TRANSCRIPT_HPP_

#include <vector>
#include <algorithm>

#include <seqan/basic.h>
#include <seqan/gff_io.h>

#define GFF_ENDPOS_SHIFT 1 // seqan stores GffRecords with a pre-shifted beginPos (0-based converted from 1-based --> if the Gff file say 8000, the GfFRecord.beginPos will say 7999)
// But the endPos are not shifted (if the Gff file says 9000, the GffRecord.endPos also says 9000)
// So, I substract 1 to endPos so that I work only with 0-based coordinates.

class exon{
public:
  exon();
  exon(int posStartInSplicedTranscript, int posStartInGenomic, int posEndInGenomic, int lg, char strand);
  void init(int posStartInSplicedTranscript, int posStartInGenomic, int posEndInGenomic, int lg, char strand);
  int getLength(void);
  int getStartCoordinateInTranscript(void);
  int getStartCoordinateGenomic(void);
  void printCoordinates(std::ostream &os);
private:
  char strand_;
  int posStartInSplicedTranscript_;
  int posStartInGenomic_;
  int posEndInGenomic_;
  int lg_;  
};


class GFFtranscriptCoordinates{
public:
  GFFtranscriptCoordinates();
  GFFtranscriptCoordinates(std::vector<seqan::GffRecord> &vr);
  bool init(std::vector<seqan::GffRecord> &vre);
  bool getAll_tToG(std::vector<int> &vpos);
  int tToG(int);
  std::string getRef(void);
  char getStrand(void);
  int getExonicTotalLength(void);
  int getStart(void);
  int getEnd(void);
  int getMinPosition(void);
  int getMaxPosition(void);
  void printStructure(std::ostream &os);
  void printSplicedSequence(std::ostream &os, seqan::Dna5String &scafSeq);
private:
  int start_;
  int end_;
  int exonsTotalLength_;
  std::string ref_; // Chromosome/scaffold name
  char strand_;
  std::vector<exon> ve_;
};



#endif
