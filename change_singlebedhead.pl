#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Script Name: change_singlebedhead.pl
# Version: 1.0
# Last Updated: 5-15-2013
#
# This script reads a bed file and changes "PercMethylation" 
# in the track name and "PercentMethylation" in the description to a new name.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input bed file
    2) NewID
" unless @ARGV == 2;

my $infile = shift(@ARGV);
my $newid = shift(@ARGV);

####################
#   Main Process   #
####################

change_singlebedhead($infile,$newid);

sub change_singlebedhead
{
	# Input
	my ($infile, $newid) = @_;
	open(IN, "<$infile") or die "Error: Cannot open $infile infile";
	# Temporary Output
	my $outfile = "temp_chr.bed";
	open(OUT, ">$outfile") or die "Error: change_bedhead.pl: cannot open $outfile outfile";
		
	# Extracting Track Name and Description
	my $firstline = <IN>;
	my @headline = split(" " , $firstline);
	my $firstcell = 0;
	while(@headline){
		# prints spaces before every cell except the first
		if ($firstcell == 0) {$firstcell = 1;} else {print OUT " ";}

		my $cell = shift(@headline);
		if($cell =~ m/name=/){print OUT "name=" , $newid;}
		elsif($cell =~ m/description=/){print OUT "description=" , $newid;}
		else {print OUT $cell}
	}
	print OUT "\n";
	
	#prints rest of file
	while(<IN>){ print OUT $_;}
	print $infile , " rename completed\n";

}

# Replace input file with output data
my $commandline = "rm " . $infile;
`$commandline`;
$commandline = "mv temp_chr.bed " . $infile;
`$commandline`;

