#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 5-28-2014
# Script Name: splitFASTAfile.pl
#
# Splits a fasta file into individual files, each with a single fasta section.
#
# Arguments: <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "splitfile.pl needs the following parameters:
    1) Input file
    2) Output files prefix (folder and prefix)
    3) Output files suffix (everything after the fasta ID, usually .fa)
" unless @ARGV == 3;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outprefix = shift(@ARGV);
my $outsuffix = shift(@ARGV);



#############
# Main Loop #
#############

while(<IN>){
	my $first = substr($_, 0, 1);
	if($first eq ">"){
		close OUT;
		my @values = split(' ', $_);
		my $FASTAname =  substr $values[0], 1;
		my $outfile = $outprefix . $FASTAname . $outsuffix;
		open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	}
	print OUT $_; 
}
close OUT;
close IN;

__END__
my $firstline = <IN>;
my @values = split(' ', $firstline);
my $FASTAname =  substr $values[0], 7;
my $outfile = $outprefix . $FASTAname . $outsuffix;
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
print OUT $firstline;

		
	print OUT $_;
	++$linecounter;
	if($linecounter > $maxlines){
		$linecounter = 1;
		close OUT;
		++$filecounter;
		if($filecounter < 10){
			$outfile = $outprefix . "00" . $filecounter . $outsuffix;
		}
		elsif($filecounter < 100){
			$outfile = $outprefix . "0" . $filecounter . $outsuffix;
		}
		else{
			$outfile = $outprefix . $filecounter . $outsuffix;
		}
		open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	}

