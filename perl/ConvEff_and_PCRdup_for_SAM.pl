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
print OUT "Sample\tConvEff\tchrM_Reads\tchrM_PCRdupsites\tchrM_PCRdupreads\tTotal_Reads\tPCRdupsites\tPCRdupreads\n";


# Global Variables
my $currentchrom = "Not Set Yet";
my %Methylation;
my %positions;
my $count = 0;


# Columns for formatting SAM files
my $chrc = 2;					#Column of chromosome in SAM file
my $startc = 3;					#Column of start location in SAM file
my $strandc = 11;				#Column of strand in SAM file



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

	my $TotalReads = 0;
	my $M_TotalReads = 0;
	my $PCRdupcount = 0;
	my $PCRdupflag = 0;
	my $PCRduploccount = 0;
	my $M_PCRdupcount = 0;
	my $M_PCRdupflag = 0;
	my $M_PCRduploccount = 0;
	
	# Previous read info
	my $prevstart = 0;
	my $prevstrand = "+";


	my $infile = shift(@ARGV);
	open(IN, "<$infile") or die "cannot open $infile infile";
	print "Starting conversion efficiency analysis of $infile \n";
	while(<IN>){
		if ($_ =~ /^@/) {next;}

		my @line = split("\t", $_);
		my $chrom = $line[$chrc];
		# Takes out weird chromosomes like chr#_random and ect.
		if ($chrom =~ /_/) {next;}

		my $start = $line[$startc];
		my $strand = substr $line[$strandc], 5,1;
		$TotalReads++;
		
		# If duplicate line
		if($prevstart == $start && $prevstrand eq $strand) {
			$PCRdupcount++;
			if($PCRdupflag == 0){
				$PCRduploccount++;
				$PCRdupflag = 1;
			}
			next;
		}
		else {$PCRdupflag = 0;}
		
		# Only look at mitochondrial chromosome
		if($chrom eq "chrM"){
			$M_TotalReads++;
			# If duplicate line
			if($prevstart == $start && $prevstrand eq $strand) {
				$M_PCRdupcount++;
				if($M_PCRdupflag == 0){
					$M_PCRduploccount++;
					$PCRdupflag = 1;
				}
			}
			else {$M_PCRdupflag = 0;}

			my $methstring = substr $line[14], 5;
			my $count = ($methstring =~ tr/x//); $CG_unmeth+=$count; $CGtot+=$count;
			$count = ($methstring =~ tr/X//); $CGtot+=$count;
			$count = ($methstring =~ tr/y//); $CHG_unmeth+=$count; $CHGtot+=$count;
			$count = ($methstring =~ tr/Y//); $CHGtot+=$count;
			$count = ($methstring =~ tr/z//); $CHH_unmeth+=$count; $CHHtot+=$count;
			$count = ($methstring =~ tr/Z//); $CHHtot+=$count;
		}
		#Update previous read information
		$prevstart = $start;
		$prevstrand = $strand;
	}
	# Print output percentages
	my $CGper = sprintf("%.4f",$CG_unmeth/$CGtot);
	my $CHGper = sprintf("%.4f",$CHG_unmeth/$CHGtot);
	my $CHHper = sprintf("%.4f",$CHH_unmeth/$CHHtot);
	my $Cper = sprintf("%.4f",($CGper+$CHGper+$CHGper) / 3);
	$PCRdupcount = $PCRdupcount + $PCRduploccount;
	$M_PCRdupcount = $M_PCRdupcount + $M_PCRduploccount;
	print OUT $infile , "\t" , $Cper ,  "\t" , $M_TotalReads , "\t" , $M_PCRduploccount , "\t" ,$M_PCRdupcount , "\t" ,$TotalReads,  "\t" , $PCRduploccount , "\t" ,$PCRdupcount , "\n";
	close IN;
}
close OUT;



