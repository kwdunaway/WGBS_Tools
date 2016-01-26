#!/usr/bin/perl 
use strict; use warnings; 
use POSIX 'floor';

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu and rgchu@ucdavis.edu
# Date: 2-24-2015
#
# This script takes windows (user defined parameters) and outputs avg meth across windows
# based on a read centric method. The script also outputs a count of CpG assays.
#
# For example: If 2 CpGs were assayed like this:
#   0.5-2
#   1-8
# The Window_permethv2 would output .75 for methylation and 2 for coverage.
# This script would output .9 for methylation and 10 for coverage.
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 usage, needs the following parameters: 
    1) Output table file
    2) Input BED file to identify what chromosomes to look at
    3) Window size
    4) Min # of CpGs per window (otherwise prints NA)
    5) Min # of reads per CpG counted
    6) Min number of files have info
    7,9+) Permeth prefix (leave off chr#.bed)
    8,10+) Name of experiments in output file
" unless @ARGV > 7;

my $outputname = shift(@ARGV);
open(OUT, ">$outputname") or die "Error: cannot open $outputname OUT outfile";
my $countoutputname = $outputname . ".count";
open(COUNT, ">$countoutputname") or die "Error: cannot open $countoutputname COUNT outfile";
my $BEDinputname = shift(@ARGV);
open(CHR, "<$BEDinputname") or die "Error: cannot open $BEDinputname bed infile";

# Size of windows in bp (ex: 20000)
my $windowsize = shift(@ARGV);
# Threshold for number of CpGs found in a window in one sample
# If it is not met, then NA is set
my $mincpg = shift(@ARGV);
# Threshold for number of reads for a CpG
# If it is not met, then NA is set
my $mincoverage = shift(@ARGV);
# Threshold for number of samples that output non-"NA" data
# If it is not met, then window is not printed to output
my $minfiles = shift(@ARGV);

# Samples
my @Permethfiles; 
my @Permethnames; 
while(@ARGV){
	push(@Permethfiles, shift(@ARGV));
	push(@Permethnames, shift(@ARGV));
}
my $commandline = "";



#######################
# Chromosome BED File #
#######################
# This is only used to identify the chromosomes
my %CHR;
while(<CHR>){
	my @line = split("\t",$_);
	if ($line[0] =~ /chr/){} else {next;}
	if ($line[0] =~ /_/){next;}	
	$CHR{$line[0]}{$line[1]} = $line[2];
}
close CHR;
print "CHRs loaded, analyzing chromosomes:\n";
foreach my $key (sort keys %CHR) {
     print "$key" , "\n";
}
print "\n";

my @samplenames;
print OUT "chr" , "\t" , "start" , "\t" , "end";
print COUNT "chr" , "\t" , "start" , "\t" , "end";
for (my $n = 0; $n < @Permethnames; $n++){
	print OUT "\t" , $Permethnames[$n];
	print COUNT "\t" , $Permethnames[$n];
}
print OUT "\n";
print COUNT "\n";



#############
# Main Loop #
#############

# For each chromosome
foreach my $key (sort keys %CHR) {
	my %MethCpG;
	my @Permeth = @Permethfiles;
	my @Names = @Permethnames;
	# For said chromosome in each sample
	while (@Permeth){
		my $inprefix = shift(@Permeth);
		my $sampname = shift(@Names);
		# Open corresponding bed file
		my $infile = $inprefix . $key . ".bed";
		# open(IN, "<$infile") or die "Error: cannot open $infile infile";
		unless(open(IN, "<$infile")) {
			$infile = $inprefix . $key . ".bed.gz";
			open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
		}
		print "Analyzing $sampname $key \n";
		while(<IN>){
			my @line = split("\t",$_);
			my @methinfo = split("-",$line[3]);
			# Check for read threshold
			if ($methinfo[1] < $mincoverage) {next;}
			my $weightedmeth = $methinfo[0] * $methinfo[1];
			my $start = $line[1];
			# Make a key for this window
			my $newkey = floor($start/$windowsize);
			# Add CpG data to hash
			# If already in hash, add to data
			if(defined $MethCpG{$sampname}{$newkey}{"sum"}) {
				$MethCpG{$sampname}{$newkey}{"sum"} = $MethCpG{$sampname}{$newkey}{"sum"} + $weightedmeth;
				$MethCpG{$sampname}{$newkey}{"coverage"} = $MethCpG{$sampname}{$newkey}{"coverage"} + $methinfo[1];
				$MethCpG{$sampname}{$newkey}{"CPGcount"}++;
				# If window meets CpG threshold,
				# increment the count for samples with data
				# (used for minimum samples threshold)
				if ($MethCpG{$sampname}{$newkey}{"CPGcount"} == $mincpg){
					if(defined $MethCpG{"count"}{$newkey}) {$MethCpG{"count"}{$newkey}++;}
					else {$MethCpG{"count"}{$newkey} = 1;}
				}
			}
			# If new to hash, add to data
			else{
				$MethCpG{$sampname}{$newkey}{"sum"} = $weightedmeth;
				$MethCpG{$sampname}{$newkey}{"coverage"} = $methinfo[1];
				$MethCpG{$sampname}{$newkey}{"CPGcount"} = 1;
				# If window meets CpG threshold,
				# increment the count for samples with data
				# (used for minimum samples threshold)
				if ($MethCpG{$sampname}{$newkey}{"CPGcount"} == $mincpg){
					if(defined $MethCpG{"count"}{$newkey}) {$MethCpG{"count"}{$newkey}++;}
					else {$MethCpG{"count"}{$newkey} = 1;}
				}
			}
		}
		close IN;		
	}	

	# Print to outfile and countfile
	@Permeth = @Permethfiles;
	@Names = @Permethnames;
	print "Printing $key \n\n";
	foreach my $countkey (sort { $a <=> $b } keys( %{$MethCpG{"count"}} ) ){
		# If window does not meet minimum sample threshold, do not print
		if ($MethCpG{"count"}{$countkey} < $minfiles) {next;}
		my $start = $countkey * $windowsize;
		my $end = $start + $windowsize;
		print OUT $key , "\t" , $start , "\t" , $end;
		print COUNT $key , "\t" , $start , "\t" , $end;
		for (my $n = 0; $n < @Names; $n++){
			# Check if window is defined
			if(defined $MethCpG{$Names[$n]}{$countkey}){
				# If # CpGs meet minimum CpG threshold, print
				if ( $MethCpG{$Names[$n]}{$countkey}{"CPGcount"} >= $mincpg) {
					my $avgmeth = sprintf ("%.3f", ($MethCpG{$Names[$n]}{$countkey}{"sum"} / $MethCpG{$Names[$n]}{$countkey}{"coverage"}));
					print OUT "\t" , $avgmeth;
					print COUNT "\t" , $MethCpG{$Names[$n]}{$countkey}{"coverage"};
				}
				# else does not meet min CpGs, print NA
				else{
					print OUT "\t" , "NA";
					print COUNT "\t" , "NA";
				}
			}
			# else print NA if no data
			else{
				print OUT "\t" , "NA";
				print COUNT "\t" , "NA";
			}
		}
		print OUT "\n";
		print COUNT "\n";
	}	
}
close OUT;
close COUNT;

