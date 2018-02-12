#!/usr/bin/perl -w
use strict;

#################################
#    Merge Interacting Reads    #
#################################

# This script was written by Anna-Lena Kranz
# This script will take two bed files and look for matching paired ends. It will then output a paired end bed file for all reads that have 
# a partner.

# Required inputs
#A-bed file 1. 
#B-bed file 1. 

#Example input for NYU HPC:
##perl /home/mrp420/yeast/scripts/merge_interacting_reads.pl MERGEREADY${TAG2}.bed MERGEREADY${TAG4}.bed ${TAG5}_MERGED.bed

#Read in my bed files
my $file1 = shift();
my $file2 = shift();
my $output = shift();

open (INPUT1, "<", $file1) || die "Couldn't open file $file1: $!.\n";
open (INPUT2, "<", $file2) || die "Couldn't open file $file2: $!.\n";
open (OUTPUT, ">", $output) || die "Couldn't open file $output: $!.\n";

my %reads;

my $count = 0;

#This looks for the pairs. It is based on using the identifier SRR* at the start of each unique ID for sequencing runs as is standard for 
# GEO. If different identifiers are used in the bed files, then modify the SRR section of the script. 

while (<INPUT1>) {
    chomp();
#    my $input2 = <INPUT2>;
#    chomp($input2);
    if (/^(SRR.*?)\t/) {
        my $header1 = $1;
        my @line1 = split("\t", $_);
        unless ($reads{$header1}) {
            $reads{$header1} = $line1[1]."\t".$line1[2]."\t".$line1[3];
        } else {
            print "again: $header1\n"; exit(0);
        }
#        my $header2;
#        if ($input2 =~ /^(SRR0381.*?)\s+/) {
#            $header2 = $1;
#        }
#        my @line1 = split("\t", $_);
#        my @line2 = split("\t", $input2);
#        unless ($header1 eq $header2) {
#            print "Problem in line $.: $line1[0] - $line2[0]!\n";
#            exit(0);
#        }
#        if ($line1[2] =~ /^chr/ && $line2[2] =~ /^chr/) {
#            print OUTPUT "$line1[0]\t$line1[2]\t$line1[3]\t$line2[0]\t$line2[2]\t$line2[3]\n";
#        }
    }   
}

#print scalar(keys(%reads))."\n";

while (<INPUT2>) {
    chomp();
    if (/^(SRR.*?)\t/) {
        my $header2 = $1;
        my @line2 = split("\t", $_);
        if ($reads{$header2}) {
            $count++;
            print OUTPUT "$header2\t".$reads{$header2}."\t$line2[1]\t$line2[2]\t$line2[3]\n";;
        }
    }  
}

#print scalar(keys(%reads))."\n";
