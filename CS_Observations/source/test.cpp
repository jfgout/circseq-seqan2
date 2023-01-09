#include "CS_ScaffoldObservations.hpp"

#include <chrono>
#include <ctime>
#include <iostream>

int main(int argc, char *argv[]){

  char tabBases[NB_BASES] = {'A', 'T', 'C', 'G'};
  
  int chrSize = 10;
  int NbObs = 80;
  std::cerr<<"Allocating the memory ...";
  CS_ScaffoldObservations sob(chrSize, "III", 2);
  std::cerr<<" done.\nPress any key to continue.";
  std::cin.get();
  std::cerr<<"\n";
  
  BasesIndices BI;
  srand(time(0));
  const clock_t begin_time = clock();
  
  for(int i=0 ; i<NbObs ; i++){

    int iStrand = rand()%2;
    int iBase = rand()%4;
    char base = tabBases[iBase];
    int pos = rand()%chrSize;

    sob.addObservation(base, 1, pos, iStrand, BI);
    
    if(i%10000000 == 0){
      std::cerr<<i<<" done ...\n";
    }
  }

  const clock_t end_time = clock();
  std::cerr<<"Time needed to add "<<NbObs<<" observations: "<<float(end_time-begin_time)/CLOCKS_PER_SEC<<"\n";

  sob.print(std::cout, true, true, true, false, BI, 0);
  std::cout<<"\n-----\nRemoving empty positions:\n";
  sob.print(std::cout, true, true, false, false, BI, 0);
  std::cout<<"\n-----\nCombining stands:\n";
  sob.print(std::cout, true, true, false, true, BI, 0);
  
  return 0;
}
