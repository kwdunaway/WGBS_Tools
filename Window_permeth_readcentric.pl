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
    2) NONE or CpG island GTF (or bed) file to mask. If no masking, put NONE
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
my $CPGinputname = shift(@ARGV);
open(CPG, "<$CPGinputname") or die "Error: cannot open $CPGinputname GTF infile";

my $windowsize = shift(@ARGV);
my $mincpg = shift(@ARGV);
my $mincoverage = shift(@ARGV);
my $minfiles = shift(@ARGV);
my @Permethfiles; 
my @Permethnames; 
while(@ARGV){
	push(@Permethfiles, shift(@ARGV));
	push(@Permethnames, shift(@ARGV));
}
my $commandline = "";



###############
# CPG Islands #
###############
# Note, at this point this is only used to identify the chromosomes
my %CPG;
while(<CPG>){
	my @line = split("\t",$_);
	if ($line[0] =~ /chr/){} else {next;}
	if ($line[0] =~ /_/){next;}	
	$CPG{$line[0]}{$line[1]} = $line[2];
}
close CPG;
print "CPGs loaded, analyzing chromosomes:\n";
foreach my $key (sort keys %CPG) {
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

foreach my $key (sort keys %CPG) {
	my %MethCpG;
	my @Permeth = @Permethfiles;
	my @Names = @Permethnames;
	while (@Permeth){
		my $inprefix = shift(@Permeth);
		my $sampname = shift(@Names);
		my $infile = $inprefix . $key . ".bed";
#		open(IN, "<$infile") or die "Error: cannot open $infile infile";
		unless(open(IN, "<$infile")) {
			$infile = $inprefix . $key . ".bed.gz";
			open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
		}
		print "Analyzing $sampname $key \n";
		while(<IN>){
			my @line = split("\t",$_);
			my @methinfo = split("-",$line[3]);
			if ($methinfo[1] < $mincoverage) {next;}
			my $weightedmeth = $methinfo[0] * $methinfo[1];
			my $start = $line[1];
			my $newkey = floor($start/$windowsize);
			if(defined $MethCpG{$sampname}{$newkey}{"sum"}) {
				$MethCpG{$sampname}{$newkey}{"sum"} = $MethCpG{$sampname}{$newkey}{"sum"} + $weightedmeth;
				$MethCpG{$sampname}{$newkey}{"coverage"} = $MethCpG{$sampname}{$newkey}{"coverage"} + $methinfo[1];
				$MethCpG{$sampname}{$newkey}{"CPGcount"}++;
				if ($MethCpG{$sampname}{$newkey}{"CPGcount"} == $mincpg){
					if(defined $MethCpG{"count"}{$newkey}) {$MethCpG{"count"}{$newkey}++;}
					else {$MethCpG{"count"}{$newkey} = 1;}
				}
			}
			else{
				$MethCpG{$sampname}{$newkey}{"sum"} = $weightedmeth;
				$MethCpG{$sampname}{$newkey}{"coverage"} = $methinfo[1];
				$MethCpG{$sampname}{$newkey}{"CPGcount"} = 1;
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
		if ($MethCpG{"count"}{$countkey} < $minfiles) {next;}
		my $start = $countkey * $windowsize;
		my $end = $start + $windowsize;
		print OUT $key , "\t" , $start , "\t" , $end;
		print COUNT $key , "\t" , $start , "\t" , $end;
		for (my $n = 0; $n < @Names; $n++){
			if(defined $MethCpG{$Names[$n]}{$countkey}){
				if ( $MethCpG{$Names[$n]}{$countkey}{"CPGcount"} >= $mincpg) {
					my $avgmeth = sprintf ("%.3f", ($MethCpG{$Names[$n]}{$countkey}{"sum"} / $MethCpG{$Names[$n]}{$countkey}{"coverage"}));
					print OUT "\t" , $avgmeth;
					print COUNT "\t" , $MethCpG{$Names[$n]}{$countkey}{"coverage"};
				}
				else{
					print OUT "\t" , "NA";
					print COUNT "\t" , "NA";
				}
			}
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

