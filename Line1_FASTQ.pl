#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Version: 6.0
# Last Updated: 1-12-2016
#
# This script looks through all raw fastq sequencing reads and finds the reads that
# have the Line1 pattern.  Then, it quantifies methylation of these sequences across 
# the four CpG sites.
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Results Table Outfile 
    2-?) Input fastq file(s) (only uses first if Filtering)
" unless @ARGV > 2;

my $results_outfile = shift(@ARGV);	
open(RESULTS, ">$results_outfile") or die "cannot open $results_outfile outfile";
my @infiles = @ARGV;

# Line1 sequence used to analyze
# TTYGTGGTGYGTYGTTTTTTAAKTYGGTT
# ctcgtggtgcgccgtttcttaagccggtctg
# CTCGTGGTGCGCCGTTTCTTAAGCCGGTC
# CTCGTGGTGCGCCGTTTCTTAA

my $inseq =         "CTCGTGGTGCGCCGTTTCTTAAGCCG";
#                    TT  T  G  TGGTG  T  G  T  T  G  TTTTTTAAGT  T  G
my @searchterms =  ("TT[CT]GTGGTG[CT]GT[CT]GTTTTTTAAGT[CT]G");
my @captureterms = ("TT([CT])GTGGTG([CT])GT([CT])GTTTTTTAAGT([CT])G");

#my %bs_seq = BS_Convert_All($inseq);
#my $printseq = Print_BS_Sequences(\%bs_seq);
#print $printseq;


#################
# In Files Loop #
#################

my %results;
my @names = ("bsfs","bsfo","bsrs","bsro");

#Original:             	CTCGTGGTGCGCCGTTTCTTAA
#BS For Same:          	TT Y  GTGGTG Y  GT Y  GTTTTTTAA
#Searchstring:         	TT[CT]GTGGTG[CT]GT[CT]GTTTTTTAA

#BS For Opp:           	CTC R  TAATAC R  CC R  TTTCTTAA
#Searchstring:         	CTC[GA]TAATAC[GA]CC[GA]TTTCTTAA

#Reverse:              	TTAAGAAACGGCGCACCACGAG
#BS Rev Same:          	TTAAGAAA Y  GG Y  GTATTA Y  GAG
#Searchstring:         	TTAAGAAA[CT]GG[CT]GTATTA[CT]GAG

#BS Rev Opp:           	TTAAAAAAC R  AC R  CACCAC R  AA
#Searchstring:         	TTAAAAAAC[GA]AC[GA]CACCAC[GA]AA

# Print table header
print RESULTS "Sample\tcount\tMeth1\tMeth2\tMeth3\tMeth4\n";

# Run process
while(@infiles){
	my $infile = shift(@infiles);

	print "\nFiltering $infile for all matched reads\n";

	my %SNPsHash = CaptureL1SNPs($infile, @captureterms);
	my $pstring = PrintL1SNPs(\%SNPsHash);
	print RESULTS $infile , $pstring;

}
close RESULTS;

###############
# Subroutines #
###############

sub CaptureL1SNPs {
	my $infile = shift;
	# if input file is zipped, gunzip and open
	if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
	else{open(IN, "<$infile") or die "cannot open $infile infile";}

	# Initialize sequence hash
	my @sequences_array = @_;
	my %SNPsequences;
	$SNPsequences{$sequences_array[0]}{"count"} = 0;
	$SNPsequences{$sequences_array[0]}{1}{"C"} = 0;
	$SNPsequences{$sequences_array[0]}{1}{"T"} = 0;
	$SNPsequences{$sequences_array[0]}{2}{"C"} = 0;
	$SNPsequences{$sequences_array[0]}{2}{"T"} = 0;
	$SNPsequences{$sequences_array[0]}{3}{"C"} = 0;
	$SNPsequences{$sequences_array[0]}{3}{"T"} = 0;
	$SNPsequences{$sequences_array[0]}{4}{"C"} = 0;
	$SNPsequences{$sequences_array[0]}{4}{"T"} = 0;
    
	while (<IN>) {
		my $ID = $_;
		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
		my $third = <IN>;
		my $quality = <IN>;
		for(my $n = 0; $n < @sequences_array; $n++){
	    	if($seq =~ /$sequences_array[$n]/){
				# Increment count and each hit
				$SNPsequences{$sequences_array[$n]}{"count"}++;
				$SNPsequences{$sequences_array[$n]}{1}{$1} = $SNPsequences{$sequences_array[$n]}{1}{$1} + 1;
				$SNPsequences{$sequences_array[$n]}{2}{$2} = $SNPsequences{$sequences_array[$n]}{2}{$2} + 1;
				$SNPsequences{$sequences_array[$n]}{3}{$3} = $SNPsequences{$sequences_array[$n]}{3}{$3} + 1;
				$SNPsequences{$sequences_array[$n]}{4}{$4} = $SNPsequences{$sequences_array[$n]}{4}{$4} + 1;
			}
		}
	}
	close IN;
	return(%SNPsequences);
}

sub PrintL1SNPs{
	my ($seq_ref) = @_;
	my %SNPsequences = %$seq_ref;
	my $printstring = "";
	foreach my $searchstring (keys %SNPsequences) {
		# Print Sample Count
		$printstring = $printstring . "\t" . $SNPsequences{$searchstring}{"count"};
		if($SNPsequences{$searchstring}{"count"} == 0){next;}
		# For each location, print percentage methylation
		for(my $SNPloc = 1; $SNPloc < 5; $SNPloc++){
			$printstring = $printstring . "\t" . sprintf("%.3f", $SNPsequences{$searchstring}{$SNPloc}{"C"}/$SNPsequences{$searchstring}{"count"});
		}
		$printstring = $printstring . "\n";
	}
	return($printstring);
}
