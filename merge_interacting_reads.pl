#!/usr/bin/perl -w
use strict;

#Input files
my $file1 = shift();
my $file2 = shift();
my $output = shift();

open (INPUT1, "<", $file1) || die "Couldn't open file $file1: $!.\n";
open (INPUT2, "<", $file2) || die "Couldn't open file $file2: $!.\n";
open (OUTPUT, ">", $output) || die "Couldn't open file $output: $!.\n";

my %reads;

my $count = 0;

#Breakdown file 1
while (<INPUT1>) {
    chomp();
    if (/^(SRR.*?)\t/) {
        my $header1 = $1;
        my @line1 = split("\t", $_);
        unless ($reads{$header1}) {
            $reads{$header1} = $line1[1]."\t".$line1[2]."\t".$line1[3];
        }
    }
}

#Breakdown file 2 and output matches
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

print $count." reads found in both pairs.\n";

close INPUT1;
close INPUT2;
close OUTPUT;
