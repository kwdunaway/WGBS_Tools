#!/usr/bin/perl 
use strict; use warnings; 

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 2-11-2016
#
# This script takes multiple Percent Methylation files and creates a single bedGraph that
# contains all of the information for those files.
#
# For example: If 2 CpGs were assayed like this:
#   chr1	2045	2046	0-2	0	+	0	0	0,0,0
#   chr1	3092	3093	1-8	0	+	0	0	210,27,27
# The script would put them into this format:
# 	chr1	2045	2046	0.5
# 	chr1	3092	3093	1
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 usage, needs the following parameters: 
    1) genome (hg38, mm10, rn6) (for chr names)
    2) Input prefix (leave off chr#.bed)
    3) Outputfile name (.bedGraph)
" unless @ARGV == 3;

my $genome = shift(@ARGV);
my $inprefix = shift(@ARGV);
my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";

my $commandline = "";



########################
# Get Chromosome Names #
########################

my %Chromosomes;
if($genome eq "hg38"){
	for(my $n = 1; $n < 23; $n++){
		my $name = "chr" . $n;
		$Chromosomes{$name} = 0;
	}
	$Chromosomes{"chrX"} = 0;
	$Chromosomes{"chrY"} = 0;
	$Chromosomes{"chrM"} = 0;
}elsif($genome eq "mm10"){
	for(my $n = 1; $n < 20; $n++){
		my $name = "chr" . $n;
		$Chromosomes{$name} = 0;
	}
	$Chromosomes{"chrX"} = 0;
	$Chromosomes{"chrY"} = 0;
	$Chromosomes{"chrM"} = 0;
}elsif($genome eq "rn6"){
	for(my $n = 1; $n < 21; $n++){
		my $name = "chr" . $n;
		$Chromosomes{$name} = 0;
	}
	$Chromosomes{"chrX"} = 0;
	$Chromosomes{"chrY"} = 0;
	$Chromosomes{"chrM"} = 0;
}
else{die "$genome is not recognized as hg38, mm10, or rn6";}


#############
# Main Loop #
#############

# For each chromosome
foreach my $cur_chr (sort keys %Chromosomes) {
	print "Starting chromosome $cur_chr" , "\t";
	# Open infile
	my $infile = $inprefix . $cur_chr . ".bed";
	unless(open(IN, "<$infile")) {
		$infile = $inprefix . $cur_chr . ".bed.gz";
		open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
	}	
	
	# Remove header if necessary
	my $firstline = <IN>; 
	my @line = split("\t",$firstline);
	if($line[0] =~ /^chr/){
		print "Header not detected in first line, converting it.\n";
		my @methinfo = split("-",$line[3]);
		print OUT $line[0] , "\t" , $line[1] , "\t", $line[2] , "\t", $methinfo[0] ,"\n";		
	}
	else {print "Header detected in first line, skipping it.\n";}


	while(<IN>){
		@line = split("\t",$_);
		my @methinfo = split("-",$line[3]);
		print OUT $line[0] , "\t" , $line[1] , "\t", $line[2] , "\t", $methinfo[0] ,"\n";		
	}
	close IN;
}
close OUT;
	
