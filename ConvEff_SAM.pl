#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Update Date: 6-5-2014
#
# Takes SAM output from BS_Seeker2 and finds the conversion efficiency by looking
# at the called methylation of mitochondria DNA (it should be completely unmethylated).
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "ConvEff_SAM.pl needs the following parameters:
    1) Output Stats
    2-?) Input SAM files
" unless @ARGV > 1;

my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
print OUT "Sample\tConvEff_CG\tConvEff_CHG\tConvEff_CHH\n";


#############
# Main Loop #
#############

# For each input file
while (@ARGV){
	my $CG_unmeth = 0;
	my $CGtot = 0;
	my $CHG_unmeth = 0;
	my $CHGtot = 0;
	my $CHH_unmeth = 0;
	my $CHHtot = 0;

	my $infile = shift(@ARGV);
	open(IN, "<$infile") or die "cannot open $infile infile";
	print "Starting conversion efficiency analysis of $infile \n";
	while(<IN>){
		my @line = split("\t", $_);
		my $chrom = $line[2];
		# Only look at mitochondrial chromosome
		if($chrom eq "chrM"){
			my $methstring = substr $line[14], 5;
			my $count = ($methstring =~ tr/x//); $CG_unmeth+=$count; $CGtot+=$count;
			$count = ($methstring =~ tr/X//); $CGtot+=$count;
			$count = ($methstring =~ tr/y//); $CHG_unmeth+=$count; $CHGtot+=$count;
			$count = ($methstring =~ tr/Y//); $CHGtot+=$count;
			$count = ($methstring =~ tr/z//); $CHH_unmeth+=$count; $CHHtot+=$count;
			$count = ($methstring =~ tr/Z//); $CHHtot+=$count;
		}
	}
	# Print output percentages
	my $CGper = sprintf("%.4f",$CG_unmeth/$CGtot);
	my $CHGper = sprintf("%.4f",$CHG_unmeth/$CHGtot);
	my $CHHper = sprintf("%.4f",$CHH_unmeth/$CHHtot);
	print OUT $infile , "\t" , $CGper , "\t" , $CHGper , "\t" , $CHHper , "\n";
	close IN;
}
close OUT;
