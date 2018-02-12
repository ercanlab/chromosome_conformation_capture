####################################
#  Fends assignment in C. elegans  #
####################################

# Will take bed file with mapped sites and assign them to fragment ends

# Required inputs
#A-Bed files with restriction enzyme sites 
#B-Paired end Bed file with Intra chromosomal interactions
#C-Paired end Bed file with Inter chromosomal interactionss

#Example input for NYU HPC:
##Rscript ~/worms/scripts/fends_elegans.R ${TAG6} ${TAG5}.bed ${TAG5}_inter.bed > outputFile${PBS_JOBID}.Rout 2>&1

#------------------------#
# Assign arguments       #
#------------------------#
print("Assign arguments...")
Sys.time()

args <- commandArgs(trailingOnly=TRUE)
args
name2<-strsplit(args[2], '.bed')[[1]][1]
directory<-paste("/scratch/mrp420/",name2,"/spike",sep="")
directory
setwd(directory)
getwd()

#-------------------------#
# Read in enzyme cut info #
#-------------------------#
print("Enzyme information...")
Sys.time()
restriction<-read.table(args[1], stringsAsFactors=F, header=F)
restriction[,3]<-NULL
restriction[which(restriction[,1]=="chrI"),]->chrI
restriction[which(restriction[,1]=="chrII"),]->chrII
restriction[which(restriction[,1]=="chrIII"),]->chrIII
restriction[which(restriction[,1]=="chrIV"),]->chrIV
restriction[which(restriction[,1]=="chrV"),]->chrV
restriction[which(restriction[,1]=="chrX"),]->chrX
restriction[which(restriction[,1]=="chrMtDNA"),]->chrM
list(chrI, chrII, chrIII, chrIV, chrV, chrX, chrM)->enzyme
key<-c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM")
names(enzyme)<-key

#-----------------------------#
# Read in mapped interactions #
#-----------------------------#
print("Mapped interactions...")
Sys.time()

#Intra-chromosomal
reads_intra<-read.table(args[2], stringsAsFactors=F, header=F)
print(reads_intra[1:10,])
#Inter-chromosomal
reads_inter<-read.table(args[3], stringsAsFactors=F, header=F)
print(reads_inter[1:10,])

#-----------------#
# Assign to FENDS #
#-----------------#
print("FEND assignment...")
Sys.time()

FENDS<-function(READ, ENZYME) {
  matrix(,1,4)->M

  if(READ[5]=="+"){
    min(which(ENZYME[[as.character(READ[1])]][,2]>as.numeric(READ[2])))->R
    ENZYME[[as.character(READ[1])]][R,1]->M[,1]
    ENZYME[[as.character(READ[1])]][R,2]->M[,2]
  }else{
        max(which(ENZYME[[as.character(READ[1])]][,2]<as.numeric(READ[2])))->R
        ENZYME[[as.character(READ[1])]][R,1]->M[,1]
        (ENZYME[[as.character(READ[1])]][R,2])+1->M[,2]
  }
  if(READ[6]=="+"){
    min(which(ENZYME[[as.character(READ[3])]][,2]>as.numeric(READ[4])))->R
    ENZYME[[as.character(READ[3])]][R,1]->M[,3]
    ENZYME[[as.character(READ[3])]][R,2]->M[,4]
  }else{
        max(which(ENZYME[[as.character(READ[3])]][,2]<as.numeric(READ[4])))->R
        ENZYME[[as.character(READ[3])]][R,1]->M[,3]
        (ENZYME[[as.character(READ[3])]][R,2])+1->M[,4]
  }
  M
}

#Apply to intra-chromsomal
print("...intra-chromosomal")
Sys.time()
apply(reads_intra, 1, FENDS, ENZYME=enzyme)->fends_intra
t(fends_intra)->fends_intra
as.data.frame(fends_intra, stringsAsFactors=F)->fends_intra
as.numeric(fends_intra[,2])->fends_intra[,2]
as.numeric(fends_intra[,4])->fends_intra[,4]
as.character(fends_intra[,3])->fends_intra[,3]
as.character(fends_intra[,3])->fends_intra[,3]
na.omit(fends_intra)->fends_intra
print(fends_intra[1:10,])

#Apply to inter-chromsomal
print("...inter-chromosomal")
Sys.time()
apply(reads_inter, 1, FENDS, ENZYME=enzyme)->fends_inter
t(fends_inter)->fends_inter
as.data.frame(fends_inter, stringsAsFactors=F)->fends_inter
as.numeric(fends_inter[,2])->fends_inter[,2]
as.numeric(fends_inter[,4])->fends_inter[,4]
as.character(fends_inter[,3])->fends_inter[,3]
as.character(fends_inter[,3])->fends_inter[,3]
na.omit(fends_inter)->fends_inter
print(fends_inter[1:10,])

#----------------------------#
# Remove adjacent FENDS      #
#----------------------------#
print("Remove adjacent FENDS...")
print("Remove adjacent FENDS...")
Sys.time()

fends_intra[which(fends_intra[,1]==fends_intra[,3] & abs(fends_intra[,2]-fends_intra[,4])>1),]->final_fends_intra

#----------------#
# Output data    #
#----------------#
print("Output data...")
Sys.time()

#Intra
filename_intra<-paste(name2,"_intra_fends",sep="")
write.table(final_fends_intra, filename_intra, col.names=F, row.names=F, quote=F, sep="\t")

#Inter
filename_inter<-paste(name2,"_inter_fends",sep="")
write.table(fends_inter, filename_inter, col.names=F, row.names=F, quote=F, sep="\t")

#ALL DONE
print("FINISHED!")
Sys.time()
