#!/usr/bin/perl 
use strict; use warnings;
use POSIX;

################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 4-62-2016
#
# This script is a modifications of AvgMeth.2col.pl where the output for each sample
# is given as three lines. First will be the # of meth reads, then # of total reads.
# The final line with be the % meth (essentially if you divide row 1 by row 2).
#
# The user can set thresholds for each read. The minimum CpG site threshold 
# will place an "NA" for the read for that experiment if the specified amount
# of CpG sites found in that read is not met. The minimum read threshold will
# ignore CpG sites with reads pertaining to that site lower than the specified
# threshold. The minimum file threshold is useful when multiple folders are
# input and requires percent methylation data (not "NA") for a read from 
# at least the specified number of folders. If the file threshold is not met,
# the bed line is not printed to the output.
#
# Note: This does NOT WORK for OVERLAPPING bed coordinates. 
#
################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: $0
    1) Output file
    2) Input BED or GTF File
    3) Minimum Read Threshold per CpG
    4,6+) Input Percent Methylation Folder Prefix (exclude \"chr\" from the path)
    5,7+) Input Sample Name (for header of output file)
" unless @ARGV > 4;

my $outputname = shift(@ARGV);	# Output with average percentage methylation per PMD
open(OUT, ">$outputname") or die "$0: Error: cannot open $outputname output file";
my $inputname = shift(@ARGV);	# Input BED or GTF File
open(BED, "<$inputname") or die "$0: Error: cannot open $inputname input BED file";
my $minreads = shift(@ARGV);	# Threshold for reads found at a CpG site
my @filenames;			# Input Percent Methylation Folder Prefix (Ex: /home/user/Permeth/Permeth_)
my @headernames;		# Input Percent Methylation Header for Column (Ex: Sample1)
while(@ARGV){
	push(@filenames,shift(@ARGV));
	push(@headernames,shift(@ARGV));
}

# Global Variables
my %bed_hash;	# Hash-hash-hash, Stores the {chromosome}{start}{end} of the BED file (columns 1, 2, 3)
my %bed_array;	# Hash-array-array, Sorted array of the hash 
		# $bed_array{chromosome}[0][2] is the line count for that chromosome
		# $bed_array{chromosome}[line count][0 for start, 1 for end]

##########################################################################
#                       Reading Input BED File                           #
# First, the BED file is loaded onto a hash by first checking for        #
# a header and then splitting every line by chromosome, start, and stop. #
# The hash is then sorted into an array to be used in later steps for    #
# faster data comparison.                                                #
##########################################################################

my $firstline = <BED>;		# Check for header
if($firstline =~ /start/){
	 print $firstline , "Header found, skipping first line!\n";
}elsif ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
	print "No header found, processing first line.\n";
	chomp($firstline);
	my @line = split("\t",$firstline);
	$bed_hash{$line[0]}{$line[1]}{$line[2]} = 1;	# Push number (no name) to hash
}
else { # If first line IS a header line
	print "Header found, skipping first line!\n";
}

# Process entire BED file
while(<BED>){ 
	chomp;
	my @line = split("\t",$_);
	if ($line[0] =~ /_/){next;}	# Ignore non-standard chromosomes
	# Check if duplicate
	if(!defined $bed_hash{$line[0]}{$line[1]}{$line[2]}) {
		$bed_hash{$line[0]}{$line[1]}{$line[2]} = 1;	# Push number (no name) to hash
	}
	# Otherwise duplicate
	else {
		print "Warning: duplicate found, skipping: $_\n";
	}
}
close BED;

# Fill in array with hash information, this is done to sort the information
foreach my $chr (keys %bed_hash){
	# Initialize line count, varies for each chromosome
	$bed_array{$chr}[0][2]=0;
	# For each position, sort numerically
	foreach my $start (sort {$a<=>$b} keys %{$bed_hash{$chr}}){
		foreach my $end (sort {$a<=>$b} keys %{$bed_hash{$chr}{$start}}){
			$bed_array{$chr}[$bed_array{$chr}[0][2]][0] = $start;	# Start
			$bed_array{$chr}[$bed_array{$chr}[0][2]][1] = $end;	# End
			$bed_array{$chr}[0][2]++;
		}
	}
}
#%bed_hash = (); # Empty the hash

print "Finished loading input BED file.\n\n";

#########################################################################
#                         Reading Input Folders                         #
# For each chromosome found in the array, search for data files for the #
# chromosome in each folder. Load the chromosome input file from each   #
# folder and fill in %outhash with information from the bed file using  #
# the %bed_array: (1) chromosome (2) start (3) end & CpG information 	#
# from the input file: (4) header name (5) percentage methylation data. #
# Print and clear %outhash for each chromosome, taking into account the #
# thresholds set by the user.                                           #
#########################################################################


# Initialize hash for storing output information
# Structure:
# $outhash{filename}{"meth"} = methylated read count
# $outhash{filename}{"total"} = total read count
my %outhash;

# Print header and initialize outhash
print OUT "Information";
for(my $i = 0; $i < @headernames; $i++){
	print OUT "\t$headernames[$i]";
	$outhash{$filenames[$i]}{"meth"} = 0;
	$outhash{$filenames[$i]}{"total"} = 0;						
}
print OUT "\n";

# Run process and print for each chromosome
foreach my $chr (sort keys %bed_array){
	print "\nLoading $chr\n";
	
	# For each sample/folder, run the chromosome bed file
	for(my $i = 0; $i < @filenames; $i++){
		my $filename = $filenames[$i] . $chr . ".bed";	# Create File Name
		my $bedpositioncounter = 0;	# Starting position for each iteration
						# Increments when previous bed array lines are no longer necessary to scan
		
		unless(open(IN, "<$filename")) {
			$filename = $filename . ".gz";
			open(IN, "gunzip -c $filename |") || die "can't open pipe to $filename";
		}
		print "Analyzing $filename\n";
		my $inline = <IN>;  	# Check for header line
		if($inline =~ /^chr/){	# Process line if data
			my @inlinearray = split("\t",$inline);
			# Check the bed array for methylation information
			for(my $k = $bedpositioncounter; $k < $bed_array{$chr}[0][2]; $k++){
				# Case 1: CpG Site is before bed array region
				# If end of in line is less than start of bed line, go to next in line
				if($inlinearray[2] < $bed_array{$chr}[$k][0]){last;}
				# Case 2: CpG Site is after bed array region
				# If end of bed line is less than start of in line, go to next bed line
				# Increment starting position, previous lines no longer necessary for further in lines
				if($bed_array{$chr}[$k][1] < $inlinearray[1]){
					$bedpositioncounter++;
					next;
				}
				# Case 3: CpG Site is within bed array region
				# If in line position lies within bed line position, add methylation information if sufficient
				if($inlinearray[1] >= $bed_array{$chr}[$k][0] && $inlinearray[2] <= $bed_array{$chr}[$k][1]){
					my @CpGmethylation = split("-", $inlinearray[3]);
					# Check if above read threshold
					my $totalreadcount = $CpGmethylation[1];
					if ($CpGmethylation[1] >= $minreads) {
						my $methreadcount = floor(($CpGmethylation[0] * $totalreadcount) + .5);
						$outhash{$filenames[$i]}{"meth"} += $methreadcount;
						$outhash{$filenames[$i]}{"total"} += $totalreadcount;						
					}
				}				
			}
		}

		while(<IN>){
			my @inlinearray = split("\t",$_);
			# Check the bed array for methylation information
			for(my $k = $bedpositioncounter; $k < $bed_array{$chr}[0][2]; $k++){
				# Case 1: CpG Site is before bed array region
				# If end of in line is less than start of bed line, go to next in line
				if($inlinearray[2] < $bed_array{$chr}[$k][0]){last;}
				# Case 2: CpG Site is after bed array region
				# If end of bed line is less than start of in line, go to next bed line
				# Increment starting position, previous lines no longer necessary for further in lines
				if($bed_array{$chr}[$k][1] < $inlinearray[1]){
					$bedpositioncounter++;
					next;
				}
				# Case 3: CpG Site is within bed array region
				# If in line position lies within bed line position, add methylation information if sufficient
				if($inlinearray[1] >= $bed_array{$chr}[$k][0] && $inlinearray[2] <= $bed_array{$chr}[$k][1]){
					my @CpGmethylation = split("-", $inlinearray[3]);
					# Check if above read threshold
					my $totalreadcount = $CpGmethylation[1];
					if ($CpGmethylation[1] >= $minreads) {
						my $methreadcount = floor(($CpGmethylation[0] * $totalreadcount) + .5);
						$outhash{$filenames[$i]}{"meth"} += $methreadcount;
						$outhash{$filenames[$i]}{"total"} += $totalreadcount;						
					}
				}
			}
		}
		close IN;
	}
}

# Print results
print "Printing results\n";
print OUT "Methylated Reads";
for(my $i = 0; $i < @filenames; $i++) {
	print OUT "\t" , $outhash{$filenames[$i]}{"meth"};
}
print OUT "\n";
print OUT "Total Reads";
for(my $i = 0; $i < @filenames; $i++) {
	print OUT "\t" , $outhash{$filenames[$i]}{"total"};
}
print OUT "\n";
print OUT "Percent Methylation";
for(my $i = 0; $i < @filenames; $i++) {
	my $average = sprintf("%.3f",$outhash{$filenames[$i]}{"meth"} / $outhash{$filenames[$i]}{"total"});
	print OUT "\t" , $average;
}
close OUT;
