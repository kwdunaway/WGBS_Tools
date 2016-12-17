#!/usr/bin/perl 
use strict; use warnings;
use POSIX 'floor';

################################################################################
# Author: Roy Chu and Keith Dunaway
# Email: rgchu@ucdavis.edu kwdunaway@ucdavis.edu
# Date: 2-24-2016
#
# This script windows sam data for coverage analysis. Note, it is important that
# coverage of your experimental samples are at least as much (if not more) than
# your control samples.
#
# Arguments:
#    <see below>
#
################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage:$0
    1) Output File
    2+) Input SAM File(s)
" unless @ARGV > 1;

my $outputfile = shift(@ARGV); # Name of output file
open(OUT, ">$outputfile") or die "I/O Error: cannot open $outputfile output file";
my @inputsam = @ARGV;

# Global variables
my %windowtable;
# Format windowtable{chr}{start}{sample} = count;
my $filecount = scalar @ARGV;
my $sizeinputsam = $#inputsam+1;
my @sampname = ();

# Header
print "Initializing.\n";
my $header = "chr";
while(@ARGV){
	my $filename = shift(@ARGV);
	$filename =~ s/.*\///;		# Remove path
	$filename =~ s/\.[^.]+$//;	# Remove extension
	$header = $header . "\t" . $filename;
	push(@sampname,$filename);
	$windowtable{"total"}{$filename} = 0;
}
$header = $header . "\n";
print OUT $header;



###############
#  Processing #
###############

for(my $n = 0; $n < @inputsam; $n++){
	my $inputfile = $inputsam[$n];
	print "Processing $inputfile\n";
	open(SAM, "<$inputfile") or die "$0 cannot open $inputfile input file";
	
	my $prevchr = "Not Set Yet";
	my $prevstart = 0;
	my $prevstrand = "+";

	while(<SAM>){
		chomp;
		my @line = split("\t",$_);
		if ($line[0] =~ /@/) {next;} #get rid of headers
		my $chrom = $line[2];
		my $start = $line[3];
		my $strand = substr $line[11], 5,1;
		if ($chrom =~ /_/) {next;} #takes out weird chromosomes like chr#_random and ect.
		if($prevchr eq $chrom && $prevstart == $start && $prevstrand eq $strand) {next;} # Skips duplicate reads (minimizes PCR bias)
		$prevchr = $chrom;
		$prevstart = $start;
		$prevstrand = "+";
		
		$windowtable{"total"}{$sampname[$n]}++;
		if(defined $windowtable{$chrom}{$sampname[$n]}){
			$windowtable{$chrom}{$sampname[$n]}++;
		}else{
			$windowtable{$chrom}{$sampname[$n]} = 1;
		}

	}
	close SAM;
}
print OUT "\n";



#######################
#  Printing OUT table #
#######################

print "Printing OUT table\n";
foreach my $chrom (sort keys %windowtable){
	print OUT $chrom;
	for(my $n = 0; $n < @sampname; $n++){
		print OUT "\t" , $windowtable{$chrom}{$sampname[$n]};
	}
	print OUT "\n";			
}
