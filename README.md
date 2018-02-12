# Chromosome Conformation Capture
A variety of methods to analyze genome-wide chromosome conformation capture techniques, such as Hi-C and 3C-seq. 


## homebrew_hic_celegans.sh
This script will take paired end data from HiC experiments and prodce a list of every pairwise intraction between restriction fragment ends (fends). This is a multi-step process:
1) Iterative mapping of each fastq file individually
2) Convert mapped files to bed format
3) Assign paired ends back together
4) Subset to inter and intra chromosomal interactions
5) Filter out interactions that are less then 1kb apart, and don't have a >1 restriction site between them. 
6) Assign each interaction to the adjacent fend
This relies on both the merge_interacting_reads.pl and fends_elegans.R script. 

## merge_interacting_reads.pl
This script was written by Anna-Lena Kranz. This script will take two bed files and look for matching paired ends. It will then output a paired end bed file for all reads that have a partner.

## fends_elegans.R
Will take bed file with mapped sites and assign them to fragment ends. This requires a bed file containing position of all the restriction sites.

## insilico_4C_celegans.R
This script will look in the two sets of chromosome conformation data produced by the Meyer lab (2015 and 2017) and export plots that show interaction counts with a specific region of interest (4C analysis). Chromosome conformation data has already been processed using the homebrew pipeline to produce a list of every interaction mapped to their fragment ends.

