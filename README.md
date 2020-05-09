# kmerCov
## What it does
The tools creates plots showing (from the top to the bottom) 

(1) the coverage of the hyper variable
genes specific kmers from the output paired end reads, 

(2) the positions, along the sequence of the specific kmers along the sequence of the detected
genotypes

(3) The proportion of the genotypes calculated as the average coverage of the kmers.

## How it works
An example command will be

python3 kmerCov.py read_1.fastq read_2.fastq condaDir cutoff

with condaDir being the full path the (mini)conda environment and cutoff is the minimum number
of reads for which a genotype is detected


