#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Last Update Date: 11-11-2014
# Version: 2.1
#
# This script trims adapter sequence from a fastq file. Currently, it only takes 
# the first 10 bases of the adapter sequence, searches for it, then trims any read
# that has a full match for the adapter.
#
# This script also "chews" back (removes) the last X (user defined) bases at the 3' 
# end. This was included because we found many bases that close to the end to have
# a unmethylated skew. 
#
# The minimum read length checks after trimming and chewing occurs. If the read does
# not meet this length, it gets removed from the dataset.
#
# There is also a filtering of reads with a quality score containing "####". This
# was done to quickly remove low quality reads that usually are a result of 
# adapter dimers.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "Usage: $0 needs the following parameters:
    1) Input FASTQ file to be trimmed (can read uncompressed or gzipped files)
    2) Output FASTQ file name (Uncompressed unless .gz is at the end of the file name)
    3) Minimum read length after adapter trimming (ex: 30)
    4) Chew back length (ex: 5 or 10)
" unless @ARGV == 4;

# Take in arguments
my $infile = shift(@ARGV);
my $outfile = shift(@ARGV);
my $minreadlength = shift(@ARGV);
my $trimlength = shift(@ARGV);

# Adapter sequence to find

#Solexa_reverse_contam
#AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAA
#my $adapter =    "AGATCGGAAG";
my $fulladapter = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAA";
print "The adapter sequence is $fulladapter\n";
print "To change, go to adapter_split.pl and change \$fulladapter\n\n";


#################
# Main Commands #
#################

# Check if zipped
my $zipflag = "N";
if ($outfile =~ /\.gz$/) {
	$outfile = substr($outfile, 0 , -3);
	$zipflag = "Y";
}
# Run process
Print_Trimmed_FASTQ_file($infile, $outfile, $fulladapter, $minreadlength, $trimlength);
# Zip again if necessary
if($zipflag eq "Y"){
	my $commandline = "gzip $outfile";
	`$commandline`;
}



###############
# Subroutines #
###############

sub Print_Trimmed_FASTQ_file{
	my ($infile, $outfile , $fulladapter, $minreadlength, $trimlength) = @_;
	
	my $adapter_seq = substr($fulladapter, 0 , 10);
	my $counter = 0;
	# Holds all sections of adapter if at the end of the sequence
	my @adapter_ends;
	my $adapter_length = length($adapter_seq);
	for(my $a = $adapter_length - 1; $a > 4; $a--) {push(@adapter_ends,substr($adapter_seq, 0, $a));}

	# File handles
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	#Opens file if gzipped
	if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";}
	#Or not gzipped
	else{open(IN, "<$infile") or die "cannot open $infile infile";}

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
		# Quality check
		if($quality =~ m/####/) {next;}
		# If adapter sequence is found
		if($seq =~ m/$adapter_seq/) {
			my @trimmedseq = split($adapter_seq, $seq);
			# Trim sequence 
			$seq = $trimmedseq[0];
			my $seqlength = length($seq);
			$quality = substr($quality, 0, length($seq));			
		}
		# Otherwise
		else{
			# Match the read ends to different lengths of the
			# adapter sequence
			for(my $t = 0; $t < @adapter_ends; $t++){
				my $adap_end_length = -1 * length($adapter_ends[$t]);
				my $end = substr($seq, $adap_end_length);
				# If match, trim
				if($end eq $adapter_ends[$t]){
					$seq = substr($seq, 0, $adap_end_length);
					$quality = substr($quality, 0, $adap_end_length);
					$t = @adapter_ends;
				}
			}
		}
		my $seqlength = length($seq);
		my $newseqlength = length($seq) - $trimlength;
		$seq = substr($seq, 0, $newseqlength);
		$quality = substr($quality, 0, $newseqlength);
		# Check if passes minimum read length
		if($seqlength < $minreadlength) {next;}
		# Print to output new sequence
		print OUT $ID , "\n" , $seq , "\n" , $third , "\n" , $quality , "\n";
	}
	close IN;
	close OUT;
	return();
}
