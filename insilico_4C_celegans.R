#################################
#    in-silico 4C C. elegans    #
#################################

# This script will look in the two sets of chromosome conformation data produced by the Meyer lab (2015 and 2017) and export plots
# that show interaction counts with a specific region of interest (4C analysis). 

# Chromosome conformation data has already been processed using the homebrew pipeline to produce a list of every interaction mapped to 
# their fragment ends.

# Required inputs
#A-Left boundary of region. 
#B-right boundary of region. 
#C-chromosome
#D-binsize
#E-Output directory
#F-Output ID
#G-Fragment

#Example input for NYU HPC:
##Rscript ~/worms/scripts/insilico_4C.R 806676 806677 chrX 5000 /scratch/mrp420/ rex40 fragment > /scratch/mrp420/reports/insilico_output.Rout 2>&1

#Arguments reported to ouptut file.
print('Script started')
Sys.time()
args <- commandArgs(trailingOnly=TRUE)
print(paste0('These are the args:',args))
query_region<-args
query_region[c(1,2,4)]<-as.numeric(query_region[c(1,2,4)])

#Test inputs - used for troubleshooting the script
#query_region<-as.data.frame(c(806676,806677,'chrX',5000,'/scratch/mrp420/insilico/','rex40'), stringsAsFactors = F)

#Create output directory
print('Make directories')
dir.create(paste0(query_region[5],'/'))
dir.create(paste0(query_region[5],query_region[6],'/'))

#Read in the data files - these are the list of fends from the meyer data. These directories are avaliable to members of the Ercan lab.  
#The equivalent files can be generated by running homebrew piepline on the meyer data.
print('Read in the data')
N2B12015<-read.table('/scratch/cgsb/ercan/GEO/2015_meyer/interactions/N2B12015_intra_fends', stringsAsFactors=F)
N2B22015<-read.table('/scratch/cgsb/ercan/GEO/2015_meyer/interactions/N2B22015_intra_fends', stringsAsFactors=F)
SDC2B12015<-read.table('/scratch/cgsb/ercan/GEO/2015_meyer/interactions/SDC2B12015_intra_fends', stringsAsFactors=F)
SDC2B22015<-read.table('/scratch/cgsb/ercan/GEO/2015_meyer/interactions/SDC2B22015_intra_fends', stringsAsFactors=F)

N2B1<-read.table('/scratch/cgsb/ercan/GEO/2017_meyer/interactions/N2B1_intra_fends', stringsAsFactors=F)
N2B2<-read.table('/scratch/cgsb/ercan/GEO/2017_meyer/interactions/N2B2_intra_fends', stringsAsFactors=F)
SDC2B1<-read.table('/scratch/cgsb/ercan/GEO/2017_meyer/interactions/SDC2B1_intra_fends', stringsAsFactors=F)
SDC2B2<-read.table('/scratch/cgsb/ercan/GEO/2017_meyer/interactions/SDC2B2_intra_fends', stringsAsFactors=F)
print('Data read-in')

#Counts for how many reads in each dataset
N2B12015_sum<-nrow(N2B12015)
N2B22015_sum<-nrow(N2B22015)
SDC2B12015_sum<-nrow(SDC2B12015)
SDC2B22015_sum<-nrow(SDC2B22015)

N2B1_sum<-nrow(N2B1)
N2B2_sum<-nrow(N2B2)
SDC2B1_sum<-nrow(SDC2B1)
SDC2B2_sum<-nrow(SDC2B2)

#Ercan lab defined Rex sites are read in for graphing purposes
rex_sites<-read.table('/scratch/cgsb/ercan/Defined_regions/rex_sites', stringsAsFactors=F)
print('Input files opened')

#Is the binsize based on a single fragment or a bin size.
print('Check what type of bait desired')
if (length(query_region)==7){
  if (query_region[7]=='fragment'){
    restriciton_sites<-read.table('/scratch/cgsb/ercan/Defined_regions/CelegansMboI.bed', stringsAsFactors=F)
    chr_specific<-restriciton_sites[which(restriciton_sites[,1]==query_region[3]),]
    index<-max(which(chr_specific[,3]<as.numeric(query_region[1]))) 
    left_boundary<-(chr_specific[index,3])-1
    right_boundary<-(chr_specific[index+1,2])+1
    
    #Define an output name
    output_name<-paste0(query_region[6],'/',query_region[6],'_fragmentbait_',query_region[4],'bins')
    
  }else{
    
    #Define the region of interest based on the binsize. This utilises the midpoint. Some modification to this definition of the bins might be 
    #required if you want to look for a large region that is bigger then a bin. 
    left_boundary<-as.numeric(query_region[1])-(as.numeric(query_region[4])/2)
    right_boundary<-as.numeric(query_region[2])+(as.numeric(query_region[4])/2)
    
    #Define an output name
    output_name<-paste0(query_region[6],'/',query_region[6],'_',query_region[4],'bins')
    
  }}

if (length(query_region)!=7){
  #Define the region of interest based on the binsize. This utilises the midpoint. Some modification to this definition of the bins might be 
  #required if you want to look for a large region that is bigger then a bin. 
  left_boundary<-as.numeric(query_region[1])-(as.numeric(query_region[4])/2)
  right_boundary<-as.numeric(query_region[2])+(as.numeric(query_region[4])/2)
  
  #Define an output name
  output_name<-paste0(query_region[6],'/',query_region[6],'_',query_region[4],'bins')
}

print(paste0('Output name is:', output_name))

#Create a function that will run on each dataset, to create a matrix of interactions. These outputs will then be used to graph my 
#viewpoint data
slid_bins<-function(input_table,input_sums) {
  
  input_name<-deparse(substitute(input_table))
  
  #Output ID
  output_ID<-paste0(query_region[6],'_',input_name,'_',query_region[4],'bins')
  
  #Find all interactions from your reigon of interest
  left_hits<-which(input_table[,1]==query_region[3]&input_table[,2]>left_boundary&input_table[,2]<right_boundary)
  right_hits<-which(input_table[,3]==query_region[3]&input_table[,4]>left_boundary&input_table[,4]<right_boundary)
  
  #Remove all in which both ends are in your region of interest
  left_partner<-input_table[left_hits,4]
  right_partner<-input_table[right_hits,2]
  
  if(length(which(left_partner>left_boundary&left_partner<right_boundary))>0){
    left_ROI<-left_partner[-which(left_partner>left_boundary&left_partner<right_boundary)]
  }else{
    left_ROI<-left_partner}
  
  if(length(which(right_partner>left_boundary&right_partner<right_boundary))>0){
    right_ROI<-right_partner[-which(right_partner>left_boundary&right_partner<right_boundary)]
  }else{
    right_ROI<-right_partner}
  
  #Add both potential directions of inteaction together.
  total_ROI<-c(left_ROI,right_ROI)
  print('Region of interest defined')
  
  #Counts are generated for each 10Kb window
  counted_ROI<-table(floor(total_ROI/as.numeric(query_region[4])))
  
  #Manipulate the rex sites so they can be graphed 
  pos<-rex_sites[1:17,2]/as.numeric(query_region[4])
  
  #how many bins are there over the x chr
  bin_number<-17720000/as.numeric(query_region[4])
  full_table<-matrix(0,bin_number,2)
  full_table[,1]<-1:bin_number
  #get the distances that have values in the table
  my_header<-as.numeric(names(counted_ROI))
  #Put the values for the table in a dat matrix for every position.
  for (i in 1:length(my_header)){
    index<-which(full_table[,1]==my_header[i])
    full_table[index,2]<-counted_ROI[i]
  }
  
  #Normalise full table to counts per million
  full_table[,2]<-(full_table[,2]/input_sums)*as.numeric(query_region[4])
  
  #Next need to create sliding average. 
  sliding_bins<-matrix(0,(bin_number-2),2)
  sliding_bins[,1]<-2:(bin_number-1)
  
  for (i in 1:nrow(sliding_bins)){
    sliding_bins[i,2]<-(sum(full_table[i,2]+full_table[i+1,2], full_table[i+2,2]))/3
  }
  assign(paste0(input_name,'sliding_bins'), sliding_bins)
  assign(paste0(input_name,'full_table'),full_table)
  
  output<-list(get(paste0(input_name,'sliding_bins')),get(paste0(input_name,'full_table')))
  return(output)
}

#Lets use the function to create the matrices of interactions across chrX
N2B1_tables<-slid_bins(N2B1,N2B1_sum)
N2B2_tables<-slid_bins(N2B2,N2B2_sum)
SDC2B1_tables<-slid_bins(SDC2B1,SDC2B1_sum)
SDC2B2_tables<-slid_bins(SDC2B2,SDC2B2_sum)

N2B12015_tables<-slid_bins(N2B12015,N2B12015_sum)
N2B22015_tables<-slid_bins(N2B22015,N2B22015_sum)
SDC2B12015_tables<-slid_bins(SDC2B12015,SDC2B12015_sum)
SDC2B22015_tables<-slid_bins(SDC2B22015,SDC2B22015_sum)

#Names for my lists of matrices
datapoints_2017<-c('N2B1_tables','N2B2_tables','SDC2B1_tables','SDC2B2_tables')
datapoints_2015<-c('N2B12015_tables','N2B22015_tables','SDC2B12015_tables','SDC2B22015_tables')


#Mainpulate the top rex sites so they can be graphed 
final_rex_sites<-c(rex_sites[1:17,5],'rex29')
pos<-c(rex_sites[1:17,2],10756306)/as.numeric(query_region[4])


#define graph parameters
if (!require("viridis")) {
  install.packages("viridis")
  library(viridis)
}
col<-viridis(4)

#Plot out the raw counts 2015
output1<-paste0(query_region[5],output_name,'_2015_raw_counts.pdf')

pdf(output1)
par(mfrow=c(2,2))
for(i in 1:4){
  my_title<-strsplit(datapoints_2015[i],'_')[[1]][1]
  plot(get(datapoints_2015[i])[[2]][,1],log10(get(datapoints_2015[i])[[2]][,2]),main=my_title,type='l', xaxt = "n", ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3.5,0))
  segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
  axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
  text(pos+40, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
  segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
  text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 
}
dev.off()

output1<-paste0(query_region[5],output_name,'_2017_raw_counts.pdf')

pdf(output1)
par(mfrow=c(2,2))
for(i in 1:4){
  my_title<-strsplit(datapoints_2017[i],'_')[[1]][1]
  plot(get(datapoints_2017[i])[[2]][,1],log10(get(datapoints_2017[i])[[2]][,2]),main=my_title,type='l', xaxt = "n", ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3.5,0))
  segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
  axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
   text(pos+40, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
  segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
  text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 
}
dev.off()


print('Raw counts plotted')


#Output the graph of sliding average.
output2<-paste0(query_region[5],output_name,'_2015_slidingwindow.pdf')

pdf(output2)
par(mfrow=c(2,2))
for(i in 1:4){
  my_title<-strsplit(datapoints_2015[i],'_')[[1]][1]
  plot(get(datapoints_2015[i])[[1]][,1],log10(get(datapoints_2015[i])[[1]][,2]),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3.5,0))
  segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
  axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
  text(pos+40, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
  segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
  text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 
}
dev.off()

output2<-paste0(query_region[5],output_name,'_2017_slidingwindow.pdf')

pdf(output2)
par(mfrow=c(2,2))
for(i in 1:4){
  my_title<-strsplit(datapoints_2017[i],'_')[[1]][1]
  plot(get(datapoints_2017[i])[[1]][,1],log10(get(datapoints_2017[i])[[1]][,2]),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3.5,0))
  segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
  axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
  text(pos+40, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
  segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
  text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 
}
dev.off()


print('Sliding window chart plotted')
Sys.time()



#What about for normalised mean datasets?

#First do mean datasets for sliding
mean_2017_sliding_wt<-(apply(cbind(get(datapoints_2017[1])[[1]][,2], get(datapoints_2017[2])[[1]][,2]), 1,mean))
mean_2017_sliding_sdc2<-(apply(cbind(get(datapoints_2017[3])[[1]][,2], get(datapoints_2017[4])[[1]][,2]), 1,mean))

mean_2015_sliding_wt<-(apply(cbind(get(datapoints_2015[1])[[1]][,2], get(datapoints_2015[2])[[1]][,2]), 1,mean))
mean_2015_sliding_sdc2<-(apply(cbind(get(datapoints_2015[3])[[1]][,2], get(datapoints_2015[4])[[1]][,2]), 1,mean))


output2<-paste0(query_region[5],output_name,'_2017_mean_norm_slidingwindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean Counts in N2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log10(mean_2017_sliding_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

my_title<-paste0('Normalised Mean Counts in SDC2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log10(mean_2017_sliding_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

dev.off()



output2<-paste0(query_region[5],output_name,'_2015_mean_norm_slidingwindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean Counts in N2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2015[1])[[1]][,1],log10(mean_2015_sliding_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

my_title<-paste0('Normalised Mean Counts in SDC2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2015[1])[[1]][,1],log10(mean_2015_sliding_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

dev.off()

#First do mean datasets for counts
mean_2017_counts_wt<-(apply(cbind(get(datapoints_2017[1])[[2]][,2], get(datapoints_2017[2])[[2]][,2]), 1,mean))
mean_2017_counts_sdc2<-(apply(cbind(get(datapoints_2017[3])[[2]][,2], get(datapoints_2017[4])[[2]][,2]), 1,mean))

mean_2015_counts_wt<-(apply(cbind(get(datapoints_2015[1])[[2]][,2], get(datapoints_2015[2])[[2]][,2]), 1,mean))
mean_2015_counts_sdc2<-(apply(cbind(get(datapoints_2015[3])[[2]][,2], get(datapoints_2015[4])[[2]][,2]), 1,mean))


output2<-paste0(query_region[5],output_name,'_2017_mean_norm_countswindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean Counts in N2 with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log10(mean_2017_counts_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

my_title<-paste0('Normalised Mean Counts in SDC2 with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log10(mean_2017_counts_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

dev.off()



output2<-paste0(query_region[5],output_name,'_2015_mean_norm_countswindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean Counts in N2 with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2015[1])[[2]][,1],log10(mean_2015_counts_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

my_title<-paste0('Normalised Mean Counts in SDC2 with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2015[1])[[2]][,1],log10(mean_2015_counts_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 

dev.off()


##Next up logratios
#First do logratios to datasets for sliding
log2_2017_sliding<-log2(mean_2017_sliding_sdc2/mean_2017_sliding_wt)
log2_2015_sliding<-log2(mean_2015_sliding_sdc2/mean_2015_sliding_wt)

output2<-paste0(query_region[5],output_name,'_log2_mean_norm_slidingwindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Log2 Normalised Mean Counts in 2017 data with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log2_2017_sliding,main=my_title, xaxt = "n",type='l', ylab='log2(sdc2/wt)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3,3))
segments(pos, 1,pos, y1 = 1.5, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, 2,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -1,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -1.5, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2,cex=0.6, labels='Bait', srt=45) 
abline(h=0, col='gray')       

my_title<-paste0('Log2 Normalised Mean Counts in 2015 data with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log2_2015_sliding,main=my_title, xaxt = "n",type='l', ylab='log2(sdc2/wt)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3,3))
segments(pos, 1,pos, y1 = 1.5, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, 2,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -1,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -1.5, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2,cex=0.6, labels='Bait', srt=45) 
abline(h=0, col='gray')       

dev.off()

#First do logratios to datasets for sliding
log2_2017_counts<-log2(mean_2017_counts_sdc2/mean_2017_counts_wt)
log2_2015_counts<-log2(mean_2015_counts_sdc2/mean_2015_counts_wt)

output2<-paste0(query_region[5],output_name,'_log2_mean_norm_countswindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Log2 Normalised Mean Counts in 2017 data with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log2_2017_counts,main=my_title, xaxt = "n",type='l', ylab='log2(sdc2/wt)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3,3))
segments(pos, 1,pos, y1 = 1.5, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, 2,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -1,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -1.5, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2,cex=0.6, labels='Bait', srt=45) 
abline(h=0, col='gray')       

my_title<-paste0('Log2 Normalised Mean Counts in 2015 data with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log2_2015_counts,main=my_title, xaxt = "n",type='l', ylab='log2(sdc2/wt)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3,3))
segments(pos, 1,pos, y1 = 1.5, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, 2,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -1,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -1.5, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2,cex=0.6, labels='Bait', srt=45) 
abline(h=0, col='gray')       
dev.off()

#First do mean datasets for sliding
mean_2017_sliding_wt<-(apply(cbind(get(datapoints_2017[1])[[1]][,2], get(datapoints_2017[2])[[1]][,2]), 1,mean))
mean_2017_sliding_sdc2<-(apply(cbind(get(datapoints_2017[3])[[1]][,2], get(datapoints_2017[4])[[1]][,2]), 1,mean))

mean_2015_sliding_wt<-(apply(cbind(get(datapoints_2015[1])[[1]][,2], get(datapoints_2015[2])[[1]][,2]), 1,mean))
mean_2015_sliding_sdc2<-(apply(cbind(get(datapoints_2015[3])[[1]][,2], get(datapoints_2015[4])[[1]][,2]), 1,mean))


output2<-paste0(query_region[5],output_name,'_2017_mean_norm_slidingwindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean Counts in N2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log10(mean_2017_sliding_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 


my_title<-paste0('Normalised Mean Counts in SDC2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log10(mean_2017_sliding_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 


dev.off()



output2<-paste0(query_region[5],output_name,'_2015_mean_norm_slidingwindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean Counts in N2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2015[1])[[1]][,1],log10(mean_2015_sliding_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 


my_title<-paste0('Normalised Mean Counts in SDC2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2015[1])[[1]][,1],log10(mean_2015_sliding_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 


dev.off()

#First do mean datasets for counts
mean_counts_wt<-(apply(cbind(mean_2017_counts_wt, mean_2015_counts_wt), 1,mean))
mean_counts_sdc2<-(apply(cbind(mean_2017_counts_sdc2, mean_2015_counts_sdc2), 1,mean))

output2<-paste0(query_region[5],output_name,'_mean_norm_countswindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean Counts in N2 with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log10(mean_counts_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 


my_title<-paste0('Normalised Mean Counts in SDC2 with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log10(mean_counts_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 


dev.off()





mean_sliding_wt<-(apply(cbind(mean_2017_sliding_wt, mean_2015_sliding_wt), 1,mean))
mean_sliding_sdc2<-(apply(cbind(mean_2017_sliding_sdc2, mean_2015_sliding_sdc2), 1,mean))

output2<-paste0(query_region[5],output_name,'_mean_norm_slidingwindow.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Normalised Mean sliding in N2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log10(mean_sliding_wt),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 


my_title<-paste0('Normalised Mean sliding in SDC2 with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[1]][,1],log10(mean_sliding_sdc2),main=my_title, xaxt = "n",type='l', ylab='log10(CPM per 10kb bin)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))),  ylim=c(-3.5,0))
segments(pos, -1,pos, y1 = -1.3, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, -0.5,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2.9,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -3.2, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -3.3,cex=0.6, labels='Bait', srt=45) 
dev.off()

#Log 2 all amalgameted data

log2_sliding<-log2(mean_sliding_sdc2/mean_sliding_wt)
log2_counts<-log2(mean_counts_sdc2/mean_counts_wt)

output2<-paste0(query_region[5],output_name,'_log2_mean_norm.pdf')
pdf(output2)
par(mfrow=c(2,1))

my_title<-paste0('Log2 Normalised Mean Counts in all data with \n counts window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log2_counts,main=my_title, xaxt = "n",type='l', ylab='log2(sdc2/wt)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4])))), ylim=c(-3,3))
segments(pos, 1,pos, y1 = 1.5, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, 2,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -1,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -1.5, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2,cex=0.6, labels='Bait', srt=45) 
abline(h=0, col='gray')       

my_title<-paste0('Log2 Normalised Mean Counts in all data with \n sliding window using ',query_region[[6]],' bait')
plot(get(datapoints_2017[1])[[2]][,1],log2_sliding,main=my_title, xaxt = "n",type='l', ylab='log2(sdc2/wt)', xlab='Chromosome position (Mb)', xlim=c(0,((18000000/as.numeric(query_region[4]))), ylim=c(-3,3))
segments(pos, 1,pos, y1 = 1.5, col=col[2], lwd=2)
axis(1,labels=seq(0,18, by=2),at=seq(0,(18000000/as.numeric(query_region[4])),by=(2000000/as.numeric(query_region[4]))))
text(pos+20, 2,cex=0.6, labels=final_rex_sites, srt=45)
segments((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -1,(left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), y1 = -1.5, col=col[4], lwd=2)
text((left_boundary+((right_boundary-left_boundary)/2))/as.numeric(query_region[4]), -2,cex=0.6, labels='Bait', srt=45) 
abline(h=0, col='gray')       
dev.off()


print('All done! Woop woop!')
