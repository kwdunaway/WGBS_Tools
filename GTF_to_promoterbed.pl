#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Last Updated: 3/25/2015
#
# This script takes regions from a GTF or bed file and 
# creates an output with the promoter region of those positions
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Input GTF or bed file name
    2) Promoter Output file name
" unless @ARGV == 2;

my $GTFfile = shift(@ARGV);
open(GTF, "<$GTFfile") or die "cannot open $GTFfile GTF or bed file";
my $promoutfile = shift(@ARGV);
open(PROM, ">$promoutfile") or die "cannot open $promoutfile promoutfile";
#my $bodyoutfile = shift(@ARGV);
#open(BODY, ">$bodyoutfile") or die "cannot open $bodyoutfile bodyoutfile";

my $promoterstart = -500;
my $promoterend = 1500;
print "The promoter start is $promoterstart from TSS\n";
print "The promoter end is +$promoterend from TSS\n";
#my $genebodystart = 3000;
#my $genebodyminsize = 5000;

###################
# Another section #
###################

while (<GTF>){
	chomp;
	my @line = split ("\t", $_);
	# if($line[2] - $line[1] < $genebodyminsize) {next;}
    
	# Strand is forward, get promoter
	my $pstart = $line[1] + $promoterstart;
	my $pend = $line[1] + $promoterend;
	# my $bstart = $line[1] + $genebodystart;
	# my $bend = $line[2];

	# Strand is reverse, get promoter
	if($line[5] eq "-"){
		$pend = $line[2] - $promoterstart;
		$pstart = $line[2] - $promoterend;
		# $bend = $line[2] - $genebodystart;
		# $bstart = $line[1];
	}
	elsif($line[5] eq "+"){}
	else{next;}
    
	# Print to output
	print PROM $line[0] ,"\t",$pstart, "\t",$pend, "\t", $line[3], "\t",$line[4], "\t",$line[5],"\n";
	# print BODY $line[0] ,"\t",$bstart, "\t",$bend, "\t", $line[3], "\t",$line[4],"\n";
}
close GTF;
close PROM;
#close BODY;
