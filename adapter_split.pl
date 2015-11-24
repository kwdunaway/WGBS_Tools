#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Last Update Date: 11-8-2014
# Version: 1.1
#
# This script separates fastq files for those with and without adapter sequence.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "Usage: $0 needs the following parameters:
    1) Input FASTQ file to be split (can read uncompressed or gzipped files)
    2) Output FASTQ file name (no adapter) (Uncompressed unless .gz is at the end of the file name)
    3) Output FASTQ file name (with adapter) (Uncompressed unless .gz is at the end of the file name)
" unless @ARGV == 3;

my $infile = shift(@ARGV);

# if file is gzipped, gunzip it first and then open
if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
else{open(IN, "<$infile") or die "cannot open $infile infile";}

# Open "no adapter" output file
my $outfile_no = shift(@ARGV);
my $outfile_no_name = $outfile_no;
if ($outfile_no_name =~ /\.gz$/) {
	$outfile_no_name = substr($outfile_no_name, 0 , -3);
}
open(OUTN, ">$outfile_no_name") or die "cannot open $outfile_no_name outfile";

# Open "with adapter" output file
my $outfile_with = shift(@ARGV);
my $outfile_with_name = $outfile_with;
if ($outfile_with_name =~ /\.gz$/) {
	$outfile_with_name = substr($outfile_with_name, 0 , -3);
}
open(OUTA, ">$outfile_with_name") or die "cannot open $outfile_with_name outfile";

# Adapter sequence to find

#Solexa_reverse_contam
#AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAA
#my $adapter = "AGATCGGAAG";
my $fulladapter = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAA";
print "The adapter sequence is $fulladapter\n";
print "To change, go to adapter_split.pl and change \$fulladapter\n\n";



#my %Sequences_After_Trimming = Print_Trimmed_FASTQ_file($infile, $outfile, $fulladapter, $minreadlength, $trimlength);
my $adapter_seq = substr($fulladapter, 0 , 10);
my $counter = 0;



####################
# Main Infile Loop #
####################

while (<IN>) {
	chomp;
	my $ID = $_;
	my $seq = <IN>;
	chop($seq); #gets rid of return character at end of sequence
	my $third = <IN>;
	$third = "+";
	my $quality = <IN>;
	chop($quality);
	# Print processing progress
	$counter++;
	if($counter % 1000000 == 0) {print "Finished procesing read:\t" , $counter , "\n";}
	# Quality Check
	if($quality =~ m/####/) {next;}
	# If adapter sequence found, place in "with adapter" output
	if($seq =~ m/$adapter_seq/) {
		print OUTA $ID , "\n" , $seq , "\n" , $third , "\n" , $quality , "\n";
	}
	# else place in "no adapter" output
	else{
		print OUTN $ID , "\n" , $seq , "\n" , $third , "\n" , $quality , "\n";
	}
}
close IN;
close OUTA;
close OUTN;


################################
# Outfile zipping if necessary #
################################

if ($outfile_no =~ /\.gz$/) {
	my $commandline = "gzip " . $outfile_no_name;
	print "Zipping $outfile_no_name\n";
	`$commandline`;
}
if ($outfile_with =~ /\.gz$/) {
	my $commandline = "gzip " . $outfile_with_name;
	print "Zipping $outfile_with_name\n";
	`$commandline`;
}
