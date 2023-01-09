##################################################################################
# This script builds all the executables for the rolling-circle analysis pipeline
#!/bin/bash

BASE_DIR=$(pwd)

# These two paths need to be modified for your local installation:
CMAKE_PREFIX_PATH="/mnt/home/software/libraries/seqan/seqan-src/util/cmake"
SEQAN_INCLUDE_PATH="/mnt/home/software/libraries/seqan/seqan-src/include"

# First part: compiling all the internal libraries that will be used by the main rolling-circle programs:
dLibs=(BasesToIndices GFF CS_Candidate CS_Consensus CS_Observations)
for dLib in ${dLibs[@]}
do
    echo -e "\n##################\nBuilding $dLib ...\n\n"
    cd $BASE_DIR/$dLib
    if [ ! -d build ]
    then
	mkdir build
    fi
    cd build
    cmake ../source/ -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH -DSEQAN_INCLUDE_PATH=$SEQAN_INCLUDE_PATH
    make
done

# Second part: compiling the programs one by one:
declare -A dProgs
dProgs["cs-findRepeats"]="CS_FindRepeats"
dProgs["refine"]="refine"
dProgs["sort-fastq-on-bam"]="sort-fastq-on-bam"
dProgs["cs-caller"]="CS_caller"

for progName in ${!dProgs[@]}
do
    progDir=${dProgs[$progName]}
    echo -e "\n##################\nBuilding $progName from $progDir...\n\n"
    cd $BASE_DIR/$progDir
    if [ ! -d build ]
    then
	mkdir build
    fi
    cd build
    cmake ../source/ -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH -DSEQAN_INCLUDE_PATH=$SEQAN_INCLUDE_PATH
    make
    cp $progName $BASE_DIR/bin/
done
