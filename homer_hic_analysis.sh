#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=homer_hic
#SBATCH --output=/scratch/mrp420/reports/slurm_homerhic_%j.out
#SBATCH --error=/scratch/mrp420/reports/slurm_homerhic_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mrp420@nyu.edu

#------------------------------------------------------------------------------#
#                                INSTRUCTIONS                                  #
#------------------------------------------------------------------------------#

#Use the script to map fastq files from Hi-C experoments.
#Once mapping is complete matrices are built and background subtraacted models built. 
#Compartment analysis is then completed.
#Lastly, inputs for the visualisation software with Juicebox
#Homer package is used to conduct this analysis

#Argument options:
# EXPID - ID for experiment that will be used in creating output files and directory
# FQ1 - Fastq file 2
# FQ2 - Fastq file 1
# DIG - Restriction site

### EXAMPLE:
# sbatch --export \
# EXPID="meyer2017_wt_rep1",\
# FQ1='/scratch/cgsb/ercan/hic/newmeyer/SRR1665088_1.fastq',\
# FQ2='/scratch/cgsb/ercan/hic/newmeyer/SRR1665082_2.fastq',\
# DIG='GATC',\
# ~/yeast/scripts/homer_hic.sh

#------------------------------------------------------------------------------#
#                                  Functions                                   #
#------------------------------------------------------------------------------#

function elapsed_time() {
    ENDTIME=$(date +%s)

    TIME=$(($ENDTIME - $1))
    if [ $TIME -lt 60 ]
    then
        echo "$TIME sec"
    elif [ $TIME -ge 60 ]  && [ $TIME -lt 3600 ]
    then
        echo "$(($TIME / 60)) min"
    else
        echo "$(($TIME / 60 / 60)) hr"
    fi
}

function check_arg() {
    if [ -z "$1" ]
    then
        echo ">>>>> Please provide values for all required arguments"
        exit 1
    fi
}

#------------------------------------------------------------------------------#
#                                  IO checks                                   #
#------------------------------------------------------------------------------#

# Check arguments
check_arg $EXPID
check_arg $FQ1
check_arg $FQ2
check_arg $DIG

# Check input files / dirs
[ -f $FQ1 ] || { echo "Could not find file: $FQ1"; exit 1; }
[ -f $FQ2 ] || { echo "Could not find file: $FQ2"; exit 1; }

#Make output directory
output_dir=/scratch/mrp420/homer_hic/$EXPID
mkdir /scratch/mrp420/homer_hic/$EXPID

#Load in packages
module load homer/intel/4.10.1
module load bowtie2/intel/2.3.2
module load samtools/intel/1.3.1
module swap r/intel/3.4.2

#------------------------------------------------------------------------------#
#                                                                              #
#                                Run pipeline                                  #
#                                                                              #
#------------------------------------------------------------------------------#

STARTTIME=$(date +%s)
echo \
ad
"------------------------------------------------------------------------------"
echo ">>>>> Started Homer HiC analysis: $EXPID"
echo \
"------------------------------------------------------------------------------"
date

#------------------------------------------------------------------------------#
#                  Align reads to reference genome with Bowtie                 #
#------------------------------------------------------------------------------#

cd $output_dir

# Trim. This cuts the reads at restriction sites. This allows better mapping. I
# If ligation junction happen within your read then the read will have multiple distinct locations it can map to, but overall no single place. 
echo 'Start Trimming' $(date +%r)
homerTools trim -3 $DIG -mis 0 -matchStart 20 -min 20 $FQ1
homerTools trim -3 $DIG -mis 0 -matchStart 20 -min 20 $FQ2
echo 'End Trimming' $(date +%r)

# Map trimmed files to the genome seperately 
echo 'Start mapping' $(date +%r)
bowtie2 --local -p 20 -x ~/worms/genomes/WS220_B2/WS220 -U ${FQ1}.trimmed > ${output_dir}/${EXPID}_1.sam
bowtie2 --local -p 20 -x ~/worms/genomes/WS220_B2/WS220 -U ${FQ2}.trimmed > ${output_dir}/${EXPID}_2.sam
echo 'End mapping' $(date +%r)

#Clean up
rm ${FQ1}.*
rm ${FQ2}.*

#------------------------------------------------------------------------------#
# Assign mapped reads to restriction enzyme fragemnt and pair with PE partner  #
#------------------------------------------------------------------------------#

# tag assignment. Initial assignment is put in subdirectory that will be home for this unfiltered tag assignment.
mkdir unfiltered

echo 'Start tag assignment' $(date +%r)
makeTagDirectory unfiltered ${output_dir}/${EXPID}_1.sam,${output_dir}/${EXPID}_2.sam
echo 'End tag assignment' $(date +%r)

# Filter tag assignment. Filters detailed:
# -removePEbg : Removes fragments likely to derive from a contiguous piece of DNA
# -restrictionSite : Removes read that are far from a restriction enzyme site
# -removeSelfLigation : Remove reads if their ends form a self ligation with adjacent restriction sites
# -removeSpikes : Remove spikes with 5x the number of reads over the average for a 10Kb region
# -tbp 1 : Only keep unique mapped partners to remove any PCR duplicates

echo 'Filter tag assignment' $(date +%r)
cp unfiltered/* .

makeTagDirectory . -update -genome ce10 -removePEbg -restrictionSite $DIG -removeSelfLigation -removeSpikes 10000 5 -tbp 1
echo 'Filter tag assignment' $(date +%r)

#------------------------------------------------------------------------------#
#    Normalise the matrix and create a model fro background interactions       #
#------------------------------------------------------------------------------#

# Normalize the marix
echo 'Create output matrix' $(date +%r)
analyzeHiC . -res 10000 -o ${EXPID}_output_matrix.txt -cpu 8
echo 'Create output matrix' $(date +%r)

echo 'Normalize matrix with balancing' $(date +%r)
analyzeHiC . -res 10000 -balance -o ${EXPID}_balanced_output_matrix.txt -cpu 8
echo 'Normalized matrix' $(date +%r)

# Background subtraction
echo 'Binning and background subtraction' $(date +%r)
analyzeHiC . -res 10000 -bgonly -cpu 8
echo 'Binning and background subtraction' $(date +%r)

# Distance normalized
echo 'Binning and background subtraction' $(date +%r)
analyzeHiC . -res 10000 -bgonly -cpu 8 -distNorm -o ${EXPID}_distnorm_output_matrix.txt
echo 'Binning and background subtraction' $(date +%r)

#------------------------------------------------------------------------------#
#            Do PCA analysis to extract compartment coordinates                #
#------------------------------------------------------------------------------#

# PCA
echo 'Start PCA' $(date +%r)
/share/apps/homer/4.10.1/intel/bin/runHiCpca.pl pcaOut . -rpath /share/apps/r/3.4.2/intel/bin/R -res 25000 -superRes 50000 -genome ce10 -cpu 10
echo 'End PCA' $(date +%r)

# PCA
echo 'Start PCA with epigentic info' $(date +%r)
/share/apps/homer/4.10.1/intel/bin/runHiCpca.pl pcaOut_emb_H3K27ac . -rpath /share/apps/r/3.4.2/intel/bin/R -res 25000 -superRes 50000 -active ~/worms/files/MACS143_e10_N2_emb_H3K27ac_LW201_LW204_LW215_ext173_159_174_peaks_final.bed -cpu 10
/share/apps/homer/4.10.1/intel/bin/runHiCpca.pl pcaOut_L3_H3K27ac . -rpath /share/apps/r/3.4.2/intel/bin/R -res 25000 -superRes 50000 -active ~/worms/files/MACS143_e10_N2_L3_H3K27ac_FE1_CTTGTA_L001_ME5_CAGATC_L005_5054_peaks_final.bed -cpu 10
echo 'End PCA with epigentic info' $(date +%r)

#------------------------------------------------------------------------------#
#                        Make Juicebox comaptible files                        #
#------------------------------------------------------------------------------#

#Juicebox
echo 'Make juicebox file' $(date +%r)
/share/apps/homer/4.10.1/intel/bin/tagDir2hicFile.pl . -juicer auto -genome ce10 -p 10 -rpath /share/apps/r/3.4.2/intel/bin/R
echo 'End juicebox file' $(date +%r)

#------------------------------------------------------------------------------#
#                        Report Significant interactions                       #
#------------------------------------------------------------------------------#

echo 'Check for signicant reads' $(date +%r)
analyzeHiC . -res 10000 -interactions ${EXPID}_significant_interactions.txt -nomatrix
echo 'Check for signicant reads' $(date +%r)


#------------------------------------------------------------------------------#
ELAPSEDTIME=$(elapsed_time $STARTTIME)
echo
echo "-----"
echo "-----"
echo "Completed pipeline in $ELAPSEDTIME"
echo \
"------------------------------------------------------------------------------"

exit 0;


