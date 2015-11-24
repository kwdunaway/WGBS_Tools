#!/usr/bin/perl 
use strict; use warnings; 
use POSIX 'floor';
##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu and rgchu@ucdavis.edu
# Date: 2-24-2015
#
# NOTE: This is a helper script, it is not meant to be run alone.
# The script that runs this is window_percentage_methylation.pl
#
# This script takes windows (user defined parameters) and outputs avg meth across windows
# based on a read centric method. The script also outputs a count of CpG assays.
#
# For example: If 2 CpGs were assayed like this:
#   0.5-2
#   1-8
# The output would be .75 for methylation and 2 for coverage.
# This script would output .9 for methylation and 10 for coverage.
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################
die "$0 usage, needs the following parameters: 
    1) Output table file
    2) Which chromosome is this? (Ex: chr1)
    3) Window size
    4) Min # of CpGs per window (otherwise prints NA)
    5) Min # of reads per CpG counted
    6) Min number of files have info
    7,9+) Permeth prefix (leave off chr#.bed)
    8,10+) Name of experiments in output file
" unless @ARGV > 7;

my $outputname = shift(@ARGV);
my $chr = shift(@ARGV);
my $chrnum = $chr;
$chr = "chr" . $chr;
$outputname = $outputname . $chr;
open(OUT, ">$outputname") or die "$0: Error: cannot open $outputname OUT file"; 
#NOTE: Appending continuously does not work (sometimes miswrites)
#open(OUT, ">>$outputname") or die "Error: cannot open $outputname OUT outfile";
#seek OUT,0,2;  # Seek the end
my $countoutputname = $outputname . ".count";
open(COUNT, ">$countoutputname") or die "$0: Error: cannot open $countoutputname OUT file"; 
#open(COUNT, ">>$countoutputname") or die "Error: cannot open $countoutputname COUNT outfile";
#seek COUNT,0,2;  # Seek the end

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




#############
# Main Loop #
#############
	my %MethCpG;
	my @Permeth = @Permethfiles;
	my @Names = @Permethnames;
	while (@Permeth){
		my $inprefix = shift(@Permeth);
		my $sampname = shift(@Names);
		my $infile = $inprefix . $chrnum . ".bed";
#		open(IN, "<$infile") or die "Error: cannot open $infile infile";
		unless(open(IN, "<$infile")) {
			$infile = $inprefix . $chrnum . ".bed.gz";
			open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
		}
		print "Analyzing $sampname $chr \n";
		while(<IN>){
			if ($_ =~ /track/) {next;}
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
	print "Printing $chr \n\n";
	foreach my $countkey (sort { $a <=> $b } keys( %{$MethCpG{"count"}} ) ){
		if ($MethCpG{"count"}{$countkey} < $minfiles) {next;}
		my $start = $countkey * $windowsize;
		my $end = $start + $windowsize;
		print OUT $chr , "\t" , $start , "\t" , $end;
		print COUNT $chr , "\t" , $start , "\t" , $end;
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

close OUT;
close COUNT;
