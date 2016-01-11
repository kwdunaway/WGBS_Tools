#!/usr/bin/perl 
use strict; use warnings; 
use POSIX 'floor';

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu and rgchu@ucdavis.edu
# Date: 1-11-2016
#
# This script takes multiple Percent Methylation files and creates a single CpG table.
#
# For example: If 2 CpGs were assayed like this:
#   0.5-2
#   1-8
# The script would put them in a table like this:
# CpG     	Sample1
# chr1_2045	1/2
# chr1_3092	8/8
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 usage, needs the following parameters: 
    1) Output table file
    2) GTF (or bed) file to determine chromosome names
    3,4+) Permeth prefix (leave off chr#.bed)
    5,6+) Name of experiments in output file
" unless @ARGV > 5;

my $outputname = shift(@ARGV);
open(OUT, ">$outputname") or die "Error: cannot open $outputname OUT outfile";
#my $countoutputname = $outputname . ".count";
#open(COUNT, ">$countoutputname") or die "Error: cannot open $countoutputname COUNT outfile";
my $CHROMOSOMESnames = shift(@ARGV);
open(CHROMOSOMES, "<$CHROMOSOMESnames") or die "Error: cannot open $CHROMOSOMESnames GTF/bed infile";

#my $windowsize = shift(@ARGV);
#my $mincpg = shift(@ARGV);
#my $mincoverage = shift(@ARGV);
#my $minfiles = shift(@ARGV);
my @Permethfiles; 
my @Permethnames; 
while(@ARGV){
	push(@Permethfiles, shift(@ARGV));
	push(@Permethnames, shift(@ARGV));
}
my $commandline = "";



########################
# Get Chromosome Names #
########################

my %Chromosomes;
while(<CHROMOSOMES>){
	my @line = split("\t",$_);
	if ($line[0] =~ /chr/){} else {next;}
	if ($line[0] =~ /_/){next;}	
	$Chromosomes{$line[0]} = 1;
}
close CHROMOSOMES;
print "CHROMOSOMES loaded, names are:\n";
foreach my $cur_chr (sort keys %Chromosomes) {
     print "$cur_chr" , "\n";
}
print "\n";



####################
# Initiate Outfile #
####################

my @samplenames;
print OUT "CpGName" , "\t" , "chr" , "\t" , "start";
for (my $n = 0; $n < @Permethnames; $n++){
	print OUT "\t" , $Permethnames[$n];
}
print OUT "\n";



#############
# Main Loop #
#############

foreach my $cur_chr (sort keys %Chromosomes) {
	my %MethCpG;
	my @Permeth = @Permethfiles;
	my @Names = @Permethnames;
	while (@Permeth){
		my $inprefix = shift(@Permeth);
		my $sampname = shift(@Names);
		my $infile = $inprefix . $cur_chr . ".bed";
		unless(open(IN, "<$infile")) {
			$infile = $inprefix . $cur_chr . ".bed.gz";
			open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
		}
		print "Analyzing $sampname $cur_chr \n";
		while(<IN>){
			my @line = split("\t",$_);
			if ($line[0] ne $cur_chr){next;}
			my $start = $line[1];
			my @methinfo = split("-",$line[3]);
			$MethCpG{$start}{$sampname}{"meth"} = floor(($methinfo[0] * $methinfo[1]) +.5);
			$MethCpG{$start}{$sampname}{"total"} = $methinfo[1];
		}
		close IN;		
	}	

	# Print to outfile
	@Names = @Permethnames;
	print "Printing $cur_chr \n\n";
	foreach my $start (sort { $a <=> $b } keys %MethCpG){
		my $cpg_name = $cur_chr . "_" . $start;
		print OUT $cpg_name , "\t" , $cur_chr , "\t" , $start;
		for (my $n = 0; $n < @Names; $n++){
			if(defined $MethCpG{$start}{$Names[$n]}{"meth"}){
					print OUT "\t" , $MethCpG{$start}{$Names[$n]}{"meth"} , "/" , $MethCpG{$start}{$Names[$n]}{"total"};
			}
			else {print OUT "\t" , "0/0";}
		}
		print OUT "\n";
	}	
}
close OUT;

