# circseq-seqan2
Rolling circle data analysis based on seqan2 library


# Dependencies:

This pipeline requires the following dependencies to be installed for compilation to work:

- seqan2
- ZLIB
- BZip2
- OpenMP
- Boost 1.67.0 (iostreams)

# Compiling:

Edit the "build-all.sh" script to specify the location of your seqan2 install and then simply run ./build-all.sh

# Usage:

## First step: identify repeats in sequencing reads with th cs-findRepeats program.
Example use:
cs-findRepeats --R1 [path to R1.fastq.gz] --R2 [path to R2.fastq.gz] --min-size 20 --max-divergence 0.1 --min-quality 30 --base-name [sampleName] -t [NB_THREADS] --gz > find.log 2> find.err

## Second step: run kallisto on the concatenated consensus sequences:
kallisto quant -i [cdna Kallisto idx] -o kallisto --pseudobam --single -l 60 -s 30 [dcsFile]

[cdna Kallisto idx] is the index file for all transcripts.
[dcsFile] is the [sampleName]-DCS.fastq.gz file produced in the previous step.

## Third step: refine the pseudo-alignment from Kallisto:
refine -k kallisto/pseudoalignments.bam -r [cdnaFile] -g [tGFF] -t [NUMBER_OF-THREADS] > refined.bam 2> refined.log
[cdnaFile] is a fasta file containing all the cDNA sequences (the one used to produce the idx used for Kallisto)
[tGFF] is a file with the GFF information for each gene in the [cdnaFile]

## 4th step: sort the refined alignment produced in the previous step:
samtools sort -o refined_sorted.bam -t MS refined.bam

## 5th step: sort the R1 and R2 reads according to the sorted bam file produced in the previous step:
sort-fastq-on-bam -b refined_sorted.bam -1 [R1.fastq.gz] -2 [R2.fastq.gz] -o R1-sorted.fastq.gz -p R2-sorted.fastq.gz

## 6th and final step: call the high-confidence mismtaches between RNA and DNA:
cs-caller -b refined_sorted.bam -p [csFile] -1 R1-sorted.fastq.gz -2 R2-sorted.fastq.gz -r [genomeFile] -g [exonsGFF] -c calls.tab.gz -o obs.tab.gz -d [minimum repeat depth - recommended = 3] -q [minimum quality stack - recommended = 100] --max-indel-size 0 > candidates.tab 2> call.log
