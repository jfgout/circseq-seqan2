#include "GFFtranscript.hpp"
#include "gff-utils.hpp"


exon::exon(){
  ;
}

exon::exon(int posStartInSplicedTranscript, int posStartInGenomic, int posEndInGenomic, int lg, char strand){
  this->init(posStartInSplicedTranscript, posStartInGenomic, posEndInGenomic, lg, strand);
}

void exon::init(int posStartInSplicedTranscript, int posStartInGenomic, int posEndInGenomic, int lg, char strand){
  this->posStartInSplicedTranscript_ = posStartInSplicedTranscript;
  this->posStartInGenomic_ = posStartInGenomic;
  this->posEndInGenomic_ = posEndInGenomic;
  this->lg_ = lg;
  this->strand_ = strand;
}


int exon::getLength(void){ return this->lg_; }

int exon::getStartCoordinateInTranscript(void){ return this->posStartInSplicedTranscript_; }

int exon::getStartCoordinateGenomic(void){ return this->posStartInGenomic_; }

void exon::printCoordinates(std::ostream &os){
  os<<this->posStartInSplicedTranscript_
    <<"\t-\t"
    << (this->posStartInSplicedTranscript_ + this->lg_ )-1
    <<"\t("
    <<this->posStartInGenomic_<<"\t-\t"<<this->posEndInGenomic_
    <<")";


  /*
  os<<this->posStartInGenomic_<<"\t-\t";
  if(this->strand_ == '='){
    os<< (this->posStartInGenomic_ - this->lg_ );
  } else {
    os<< (this->posStartInGenomic_ + this->lg_ );
  }
  os<<")";
  */
}


////////////////////////////////////////////////////////////////////////////////


GFFtranscriptCoordinates::GFFtranscriptCoordinates(){
  ;
}

GFFtranscriptCoordinates::GFFtranscriptCoordinates(std::vector<seqan::GffRecord> &vr){
  this->init(vr);
}


bool GFFtranscriptCoordinates::init(std::vector<seqan::GffRecord> &vre){

  int nbe = vre.size();
  if( nbe < 1 ){ return false; }
  
  char strand = vre[0].strand;
  this->strand_ = strand;
  std::string ref(seqan::toCString(vre[0].ref));
  this->ref_ = ref;
  
  if( strand == '-' ){
    std::sort(vre.begin(), vre.end(), sortGffRecordDecreasingGenomic);
    this->start_ = std::max(vre[0].beginPos , vre[0].endPos-GFF_ENDPOS_SHIFT);
    this->end_ = std::min(vre[nbe-1].beginPos, vre[nbe-1].endPos-GFF_ENDPOS_SHIFT);
  } else {
    std::sort(vre.begin(), vre.end(), sortGffRecordIncreasingGenomic);
    this->start_ = std::min(vre[0].beginPos , vre[0].endPos-GFF_ENDPOS_SHIFT);
    this->end_ = std::max(vre[nbe-1].beginPos, vre[nbe-1].endPos-GFF_ENDPOS_SHIFT);
  }

  int nbr = vre.size();
  if(nbr < 1){ return false; }
  this->ve_.resize(nbr);

  

  //int minPosGenomic = std::max(vre[0].beginPos, vre[0].endPos); // Just in case I deal with a strange GFF fil that has end before start for features on the minus strand...
  int currentExonStartInTranscript = 0;
  int previousExonSize = 0;
  int totalExonsSizes = 0;
  
  for(int ie = 0 ; ie < nbr ; ie++){
    //std::cout<<"Record: "<<vre[ie].beginPos<<" - "<<vre[ie].endPos<<"\n";
    currentExonStartInTranscript += previousExonSize;
    int currentExonMinGenomicPos = std::min(vre[ie].beginPos, vre[ie].endPos-GFF_ENDPOS_SHIFT);
    int currentExonMaxGenomicPos = std::max(vre[ie].beginPos, vre[ie].endPos-GFF_ENDPOS_SHIFT);
    int currentExonSize = (currentExonMaxGenomicPos - currentExonMinGenomicPos) + 1;
    totalExonsSizes += currentExonSize;
    int currentExonStartGenomic = currentExonMinGenomicPos;
    int currentExonEndGenomic = currentExonMaxGenomicPos;
    if( strand == '-' ){
      currentExonStartGenomic = currentExonMaxGenomicPos;
      currentExonEndGenomic = currentExonMinGenomicPos;
    }

    exon ve(currentExonStartInTranscript, currentExonStartGenomic, currentExonEndGenomic, currentExonSize, strand);
    this->ve_[ie] = ve;
    previousExonSize = currentExonSize;    
  }
  this->exonsTotalLength_ = totalExonsSizes;
  return true;
}

std::string GFFtranscriptCoordinates::getRef(void){
  return this->ref_;
}

char GFFtranscriptCoordinates::getStrand(void){
  return this->strand_;
}

int GFFtranscriptCoordinates::getStart(void){
  return this->start_;
}

int GFFtranscriptCoordinates::getEnd(void){
  return this->end_;
}

int GFFtranscriptCoordinates::getMinPosition(void){
  if( this->strand_ == '-' ){ return this->end_; }
  return this->start_;
}

int GFFtranscriptCoordinates::getMaxPosition(void){
  if( this->strand_ == '-' ){ return this->start_; }
  return this->end_;
}


int GFFtranscriptCoordinates::tToG(int pos){

  if( pos<0 || pos>=this->exonsTotalLength_ ){
    return -1;
  }

  int ie = 0;
  int nbe = this->ve_.size();
  int exonStartGenomic = this->start_;
  int exonEndTranscript = this->ve_[ie].getLength() - 1;

  while( exonEndTranscript < pos && ie < (nbe-1) ){
    ie++;
    exonEndTranscript += this->ve_[ie].getLength();
    exonStartGenomic = this->ve_[ie].getStartCoordinateGenomic();
  }
  
  if( ie == nbe ){
    std::cerr<<"PROBLEM WITH MAPPING TRANSCRIPTOMIC COORDINATE To GENOMIC!\n";
    return -1;
  }
  
  int posInCurrentExon = pos - this->ve_[ie].getStartCoordinateInTranscript();
  int genomicPos = -1;
  if( this->strand_ == '-' ){
    genomicPos = this->ve_[ie].getStartCoordinateGenomic() - posInCurrentExon;
  } else {
    genomicPos = this->ve_[ie].getStartCoordinateGenomic() + posInCurrentExon;
  }
  return genomicPos;
}


bool GFFtranscriptCoordinates::getAll_tToG(std::vector<int> &vpos){
  vpos.resize(this->exonsTotalLength_);
  for(int iPos=0 ; iPos<this->exonsTotalLength_ ; iPos++){
    int gPos = this->tToG(iPos);
    if( gPos < 0 ){ return false; }
    vpos[iPos] = gPos;
  }
  return true;
}

int GFFtranscriptCoordinates::getExonicTotalLength(void){
  return this->exonsTotalLength_;
}

void GFFtranscriptCoordinates::printStructure(std::ostream &os){
  int nbe = this->ve_.size();
  for(int ie=0 ; ie<nbe ; ie++){
    this->ve_[ie].printCoordinates(os);
    os<<"\n";
  }
  
}
