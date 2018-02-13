#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=homebrew_hic
#SBATCH --output=/scratch/mrp420/reports/slurm_homebrew_%j.out
#SBATCH --error=/scratch/mrp420/reports/slurm_homebrew_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mrp420@nyu.edu

#################################
#    Homebrew HiC C. elegans    #
#################################

# This script will take paired end data from HiC experiments and prodce a list of every pairwise intraction between restriction fragment 
# ends (fends). This is a multi-step process:
# 1) Iterative mapping of each fastq file individually
# 2) Convert mapped files to bed format
# 3) Assign paired ends back together
# 4) Subset to inter and intra chromosomal interactions
# 5) Filter out interactions that are less then 1kb apart, and don't have a >1 restriction site between them. 
# 6) Assign each interaction to the adjacent fend

# Required inputs
#A-fastq file 1. 
#B-fastq file 2. 
#C-bed file containing position of all the restriction sites

#Example input for NYU HPC:
##sbatch --export arg1='/scratch/cgsb/ercan/GEO/2015_meyer/SDC2B12015_X_X_R1_X.fastq',arg2='/scratch/cgsb/ercan/GEO/2015_meyer/SDC2B12015_X_X_R2_X.fastq',arg3='/home/mrp420/worms/restrict/CelegansMboI.bed' ~/worms/scripts/homebrew_hic_elegeans_geo.sh

module load bowtie/gnu/1.2.0

#Start time
echo 'Start Time'
echo $(date +%x_%r)

# Lets get the input variables and generate new variables from them. This is specific to the naming format of fastqs that NYU gencore used 
# to give out their sequnceing runs (SEQID_X_X_R1orR2_X.fastq). Fastqs must be converted to this format, or this section can be modified 
# to accept differnt input types.

TAG1=$arg1
TAG7=$( echo ${TAG1}|awk -F '[/_]' '{print $8"_"$11}' )
TAG3=$arg2
TAG9=$( echo ${TAG3}|awk -F '[/_]' '{print $8"_"$11}' )
TAG5=$( echo ${TAG3}|awk -F '[/_]' '{print $8}' )
TAG2=$(echo ${TAG7})
TAG4=$(echo ${TAG9})

#Lets create the output directories
#This is coded to output into my directory. Needs to be recoded to output into desired location throughout the pipeline. 
cd /scratch/mrp420/homebrew_hic
mkdir ${TAG5}
mkdir ${TAG5}/${TAG2}
mkdir ${TAG5}/${TAG2}/mapping
mkdir ${TAG5}/${TAG4}
mkdir ${TAG5}/${TAG4}/mapping
TAG6=$arg3

#Output the library identifier to the error report
echo $TAG5

## 1) Iterative mapping of each fastq file individually.
# Location of the bowtie index is hardcoded here. Needs to replaced to run by others. 

#How long are the reads? I need this to be able to trim the reads down to do iterative mapping. As interactions could be present within the read, mapping the whole read could be problematic as some reads could be from several places in the genome.
TAG10=$(head -1 ${TAG1} | cut -d"=" -f2)


#Iterative mapping of the first fastq file
bowtie -q -5 1 -3 expr ${TAG10} - 50 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP50.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG1} ${TAG5}/${TAG2}/mapping/L_50.sam
bowtie -q -5 1 -3 expr ${TAG10} - 45 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP45.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG2}/mapping/L_UNMAP50.fastq ${TAG5}/${TAG2}/mapping/L_45.sam
bowtie -q -5 1 -3 expr ${TAG10} - 40 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP40.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG2}/mapping/L_UNMAP45.fastq ${TAG5}/${TAG2}/mapping/L_40.sam
bowtie -q -5 1 -3 expr ${TAG10} - 35 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP35.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG2}/mapping/L_UNMAP40.fastq ${TAG5}/${TAG2}/mapping/L_35.sam
bowtie -q -5 1 -3 expr ${TAG10} - 30 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP30.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG2}/mapping/L_UNMAP35.fastq ${TAG5}/${TAG2}/mapping/L_30.sam
bowtie -q -5 1 -3 expr ${TAG10} - 25 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP25.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG2}/mapping/L_UNMAP30.fastq ${TAG5}/${TAG2}/mapping/L_25.sam
bowtie -q -5 1 -3 expr ${TAG10} - 20 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP20.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG2}/mapping/L_UNMAP25.fastq ${TAG5}/${TAG2}/mapping/L_20.sam
bowtie -q -5 1 -3 expr ${TAG10} - 15 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG2}/mapping/L_UNMAP15.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG2}/mapping/L_UNMAP20.fastq ${TAG5}/${TAG2}/mapping/L_15.sam

echo 'worm 1 mapped'

#Iterative mapping of the second fastq file
bowtie -q -5 1 -3 expr ${TAG10} - 50-v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP50.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG3} ${TAG5}/${TAG4}/mapping/R_50.sam
bowtie -q -5 1 -3 expr ${TAG10} - 45 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP45.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG4}/mapping/R_UNMAP50.fastq ${TAG5}/${TAG4}/mapping/R_45.sam
bowtie -q -5 1 -3 expr ${TAG10} - 40 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP40.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG4}/mapping/R_UNMAP45.fastq ${TAG5}/${TAG4}/mapping/R_40.sam
bowtie -q -5 1 -3 expr ${TAG10} - 35 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP35.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG4}/mapping/R_UNMAP40.fastq ${TAG5}/${TAG4}/mapping/R_35.sam
bowtie -q -5 1 -3 expr ${TAG10} - 30 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP30.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG4}/mapping/R_UNMAP35.fastq ${TAG5}/${TAG4}/mapping/R_30.sam
bowtie -q -5 1 -3 expr ${TAG10} - 25 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP25.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG4}/mapping/R_UNMAP30.fastq ${TAG5}/${TAG4}/mapping/R_25.sam
bowtie -q -5 1 -3 expr ${TAG10} - 20 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP20.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG4}/mapping/R_UNMAP25.fastq ${TAG5}/${TAG4}/mapping/R_20.sam
bowtie -q -5 1 -3 expr ${TAG10} - 15 -v 0 -m 1 -p 8 --un ${TAG5}/${TAG4}/mapping/R_UNMAP15.fastq --seed=123 -S /home/mrp420/worms/genomes/WS220/WS220_ucsc ${TAG5}/${TAG4}/mapping/R_UNMAP20.fastq ${TAG5}/${TAG4}/mapping/R_15.sam

echo 'worm 2 mapped'

## 2) Convert mapped files to bed format

module load samtools/intel/1.3.1
module load bedtools/intel/2.26.0 

cd /scratch/mrp420/homebrew_hic/${TAG5}/${TAG2}/mapping

rm *.fastq

#Convert my files to a bed file
samtools view -b -S L_50.sam > L_50.bam
samtools view -h -F 4 -b L_50.bam > L_50_mapped.bam
bedtools bamtobed -i L_50_mapped.bam > L_50_mapped.bed

samtools view -b -S L_45.sam > L_45.bam
samtools view -h -F 4 -b L_45.bam > L_45_mapped.bam
bedtools bamtobed -i L_45_mapped.bam > L_45_mapped.bed

samtools view -b -S L_40.sam > L_40.bam
samtools view -h -F 4 -b L_40.bam > L_40_mapped.bam
bedtools bamtobed -i L_40_mapped.bam > L_40_mapped.bed

samtools view -b -S L_35.sam > L_35.bam
samtools view -h -F 4 -b L_35.bam > L_35_mapped.bam
bedtools bamtobed -i L_35_mapped.bam > L_35_mapped.bed

samtools view -b -S L_30.sam > L_30.bam
samtools view -h -F 4 -b L_30.bam > L_30_mapped.bam
bedtools bamtobed -i L_30_mapped.bam > L_30_mapped.bed

samtools view -b -S L_25.sam > L_25.bam
samtools view -h -F 4 -b L_25.bam > L_25_mapped.bam
bedtools bamtobed -i L_25_mapped.bam > L_25_mapped.bed

samtools view -b -S L_20.sam > L_20.bam
samtools view -h -F 4 -b L_20.bam > L_20_mapped.bam
bedtools bamtobed -i L_20_mapped.bam > L_20_mapped.bed

samtools view -b -S L_15.sam > L_15.bam
samtools view -h -F 4 -b L_15.bam > L_15_mapped.bam
bedtools bamtobed -i L_15_mapped.bam > L_15_mapped.bed

echo 'worm 1 file is converted'

rm *.bam
rm *.sam

#Merge the bed files
cat L_15_mapped.bed L_20_mapped.bed L_25_mapped.bed L_30_mapped.bed L_35_mapped.bed L_40_mapped.bed L_45_mapped.bed L_50_mapped.bed > ../..//${TAG2}.bed

echo 'worm 1 files are collated'

#Report out the mapping efficiencyl
rm ../../${TAG5}_report 

echo "Total mapped reads" ${TAG7} >> ../../${TAG5}_report
wc -l ../../${TAG7}.bed >> ../../${TAG5}_report

echo "Mapped reads - " >> ../${TAG5}_report
wc -l L_50_mapped.bed >> ../../${TAG5}_report  
wc -l L_45_mapped.bed >> ../../${TAG5}_report  
wc -l L_40_mapped.bed >> ../../${TAG5}_report  
wc -l L_35_mapped.bed >> ../../${TAG5}_report  
wc -l L_30_mapped.bed >> ../../${TAG5}_report  
wc -l L_25_mapped.bed >> ../../${TAG5}_report  
wc -l L_20_mapped.bed >> ../../${TAG5}_report  
wc -l L_15_mapped.bed >> ../../${TAG5}_report    

echo "Total mapped reads - " ${TAG2} >> ../..//${TAG5}_report
wc -l ../../${TAG2}.bed >> ../..//${TAG5}_report

rm *.bed

echo 'worm 1 files are counted'

cd /scratch/mrp420/homebrew_hic/${TAG5}/${TAG4}/mapping

rm *.fastq

samtools view -b -S R_50.sam > R_50.bam
samtools view -h -F 4 -b R_50.bam > R_50_mapped.bam
bedtools bamtobed -i R_50_mapped.bam > R_50_mapped.bed

samtools view -b -S R_45.sam > R_45.bam
samtools view -h -F 4 -b R_45.bam > R_45_mapped.bam
bedtools bamtobed -i R_45_mapped.bam > R_45_mapped.bed

samtools view -b -S R_40.sam > R_40.bam
samtools view -h -F 4 -b R_40.bam > R_40_mapped.bam
bedtools bamtobed -i R_40_mapped.bam > R_40_mapped.bed

samtools view -b -S R_35.sam > R_35.bam
samtools view -h -F 4 -b R_35.bam > R_35_mapped.bam
bedtools bamtobed -i R_35_mapped.bam > R_35_mapped.bed

samtools view -b -S R_30.sam > R_30.bam
samtools view -h -F 4 -b R_30.bam > R_30_mapped.bam
bedtools bamtobed -i R_30_mapped.bam > R_30_mapped.bed

samtools view -b -S R_25.sam > R_25.bam
samtools view -h -F 4 -b R_25.bam > R_25_mapped.bam
bedtools bamtobed -i R_25_mapped.bam > R_25_mapped.bed

samtools view -b -S R_20.sam > R_20.bam
samtools view -h -F 4 -b R_20.bam > R_20_mapped.bam
bedtools bamtobed -i R_20_mapped.bam > R_20_mapped.bed

samtools view -b -S R_15.sam > R_15.bam
samtools view -h -F 4 -b R_15.bam > R_15_mapped.bam
bedtools bamtobed -i R_15_mapped.bam > R_15_mapped.bed

echo 'worm 2 file is converted'

rm *.bam
rm *.sam

cat R_15_mapped.bed R_20_mapped.bed R_25_mapped.bed R_30_mapped.bed R_35_mapped.bed R_40_mapped.bed R_45_mapped.bed R_50_mapped.bed > ../..//${TAG4}.bed

echo 'worm 2 files are collated'

echo "Total mapped reads" ${TAG9} >> ../../${TAG5}_report
wc -l ../../${TAG9}.bed >> ../../${TAG5}_report

echo "Mapped reads - " >> ../${TAG5}_report
wc -l R_50_mapped.bed >> ../../${TAG5}_report  
wc -l R_45_mapped.bed >> ../../${TAG5}_report  
wc -l R_40_mapped.bed >> ../../${TAG5}_report  
wc -l R_35_mapped.bed >> ../../${TAG5}_report  
wc -l R_30_mapped.bed >> ../../${TAG5}_report  
wc -l R_25_mapped.bed >> ../../${TAG5}_report  
wc -l R_20_mapped.bed >> ../../${TAG5}_report  
wc -l R_15_mapped.bed >> ../../${TAG5}_report    

echo "Total mapped reads - " ${TAG4} >> ../../${TAG5}_report
wc -l ../../${TAG4}.bed >> ../../${TAG5}_report

rm *.bed

echo 'worm 2 files are counted'


## 3) Assign paired ends back together

cd /scratch/mrp420/homebrew_hic/${TAG5}/

awk '{OFS="\t"; print $4,$1,$2,$6}' ${TAG2}.bed > MERGEREADY${TAG2}.bed
awk '{OFS="\t"; print $4,$1,$2,$6}' ${TAG4}.bed > MERGEREADY${TAG4}.bed

module load perl/intel/5.24.0

perl /home/mrp420/yeast/scripts/merge_interacting_reads.pl MERGEREADY${TAG2}.bed MERGEREADY${TAG4}.bed ${TAG5}_MERGED.bed

echo 'paired end files are matched together'

echo "Mapped reads in both pairs - no rDNA" ${TAG5} >> ${TAG5}_report
wc -l ${TAG5}_MERGED.bed >> ${TAG5}_report

## 4) Subset to inter and intra chromosomal interactions

awk '{OFS="\t"} {if ($2 == $5) print $2,$3,$6,$4,$7,$1}' ${TAG5}_MERGED.bed > ${TAG5}_intra.bed

awk '{OFS="\t"} {if ($2 < $3) print $1,$2,$3,$4,$5,$6;
else
print $1,$3,$2,$4,$5,$6;
}' ${TAG5}_intra.bed > ${TAG5}_intrasorted.bed
awk '{OFS="\t"} {if ($2 != $5) print $2,$3,$5,$6,$4,$7,$1}' ${TAG5}_MERGED.bed > ${TAG5}_inter.bed

echo 'intra and inter interactions are seperated from each other'

echo "Intrachromsomal(Prefilter)" ${TAG5} >> ${TAG5}_report
wc -l ${TAG5}_intrasorted.bed >> ${TAG5}_report

echo "Interchromosomal" ${TAG5} >> ${TAG5}_report
wc -l ${TAG5}_inter.bed >> ${TAG5}_report

## 5) Filter out interactions that are less then 1kb apart, and don't have a >1 restriction site between them. 

bedtools intersect -u -a ${TAG5}_intrasorted.bed -b ${TAG6} > ${TAG5}_withrestrict.bed

echo 'interactions that have at least one restriction site between them are filtered out'

echo "Intrachromsomal with Restriction Site" ${TAG5} >> ${TAG5}_report
wc -l ${TAG5}_withrestrict.bed >> ${TAG5}_report

echo ${TAG6} >> ${TAG5}_report

echo "Intrachromsomal with over 1000bp" ${TAG5} >> ${TAG5}_report
awk '{OFS="\t"} {if ($3-$2 >= 1000) print $1,$2,$1,$3,$4,$5,$6}' ${TAG5}_withrestrict.bed > ${TAG5}.bed
wc -l ${TAG5}.bed >> ${TAG5}_report

echo 'interactions must be at least 1kb apart'

## 6) Assign each interaction to the adjacent fend

module load r/intel/3.3.2

Rscript ~/worms/scripts/fends_elegans.R ${TAG6} ${TAG5}.bed ${TAG5}_inter.bed > outputFile${PBS_JOBID}.Rout 2>&1

echo 'interactions are assigned to fends'

#End time
echo 'End Time'
echo $(date +%x_%r)

exit 0;
