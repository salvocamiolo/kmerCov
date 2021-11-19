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

python3 kmerCov.py read_1.fastq read_2.fastq condaDir cutoff outputFolder

with condaDir being the full path the (mini)conda environment and cutoff is the minimum number
of reads for which a genotype is detected. Two software must be installed in the conda environment:

(a) jellyfish 

conda install -c conda-forge jellyfish

or

conda install -c conda-forge/label/gcc7 jellyfish

or

conda install -c conda-forge/label/cf201901 jellyfish

or

conda install -c conda-forge/label/cf202003 jellyfish


(b) mafft

conda install -c bioconda mafft

or

conda install -c bioconda/label/cf201901 mafft




