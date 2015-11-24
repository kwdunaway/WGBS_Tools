#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 12/4/2014
#
# Takes a BSSeeker2 log file and and filters it for only the information I care about.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Outfile (txt tab delimited)
    2+) Infile(s) of logs
" unless @ARGV > 1;

my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
print OUT "FASTQ Files		 Total bases of uniquely mapped reads 		Mappability		mCG		mCHG		mCHH	
noadap	trimmed	 noadap 	 trimmed 	noadap	trimmed	noadap	trimmed	noadap	trimmed	noadap	trimmed
";
my @Infiles = @ARGV;
my @usefulinfo;

##############
# Main loops #
##############

while(@Infiles){
	my $infile = shift(@Infiles);
	open(IN, "<$infile") or die "cannot open $infile infile";
	while (<IN>)
	{
	    chomp;
	    my @line = split ("\ ", $_);
	    if($_ =~ /Read filename/){
	    	push(@usefulinfo,$line[4]);
	    }elsif($_ =~ /Mappability/){
	    	push(@usefulinfo,$line[3]);	    
	    }elsif($_ =~ /Total bases of uniquely mapped reads/){
	    	push(@usefulinfo,$line[8]);	    
	    }elsif($_ =~ /  mCG/){
	    	push(@usefulinfo,$line[3]);	    
	    }elsif($_ =~ /  mCHG/){
	    	push(@usefulinfo,$line[3]);	    
	    }elsif($_ =~ /  mCHH/){
	    	push(@usefulinfo,$line[3]);	    
	    }
	    if(@usefulinfo == 12){
	    	print OUT $usefulinfo[0] , "\t" , $usefulinfo[6] , "\t" , 
	    	$usefulinfo[2] , "\t" , $usefulinfo[8] , "\t" , 
	    	$usefulinfo[1] , "\t" , $usefulinfo[7] , "\t" , 
	    	$usefulinfo[3] , "\t" , $usefulinfo[9] , "\t" , 
	    	$usefulinfo[4] , "\t" , $usefulinfo[10] , "\t" , 
	    	$usefulinfo[5] , "\t" , $usefulinfo[11] , "\n";
	    	@usefulinfo = ();
	    }
	}
	close IN;
}
close OUT;
