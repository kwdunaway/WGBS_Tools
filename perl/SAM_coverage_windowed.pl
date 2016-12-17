#!/usr/bin/perl 
use strict; use warnings;
use POSIX 'floor';

################################################################################
# Author: Roy Chu and Keith Dunaway
# Email: rgchu@ucdavis.edu kwdunaway@ucdavis.edu
# Date: 4-8-2015
#
# This script windows sam data for coverage analysis. Note, it is important that
# coverage of your experimental samples are at least as much (if not more) than
# your control samples.
#
# Normalizing factor = Average CTRLs and divide by coverage and divide by 2
# Divide other samples by normalizing factor
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
    2) Window Size
    3) Minimum control coverage (ex: 30)
    4) Number of control samples (input those sam files first) (ex: 6)
    5+) Input SAM File(s)
" unless @ARGV > 5;

my $outputfile = shift(@ARGV); # Name of output file
open(OUT, ">$outputfile") or die "Window_Coverage_Normalized.pl Error: cannot open $outputfile output file";
my $windowsize = shift(@ARGV);
my $mincoverage = shift(@ARGV);
my $ctrlnum = shift(@ARGV);
my @inputsam = @ARGV;

# Global variables
my %windowtable;
# Format windowtable{chr}{start}{sample} = count;
my $filecount = scalar @ARGV;
my $sizeinputsam = $#inputsam+1;
my @sampname = ();

# Header
print "Initializing.\n";
my $header = "chr\tstart\tend";
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
	open(SAM, "<$inputfile") or die "Window_Coverage_Normalized.pl cannot open $inputfile input file";
	
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

		my $key = floor($start/$windowsize);
		
		if( defined $windowtable{$chrom}{$key}{$sampname[$n]}) {
			$windowtable{$chrom}{$key}{$sampname[$n]}++;
			$windowtable{"total"}{$sampname[$n]}++;
		}
		else{
			$windowtable{$chrom}{$key}{$sampname[$n]} = 1;
			$windowtable{"total"}{$sampname[$n]}++;
		}
	}
	close SAM;
}

#Calculate conversion rates
print "Calculating converian rates of all samples\n";
#my %convertedtable;
foreach my $chrom (sort keys %windowtable){
	if($chrom eq "total") {next;}
	foreach my $key (sort { $a <=> $b } keys( %{$windowtable{$chrom}} ) ){
		my $countflag = 0;
		my $normfactor = 0;
		for(my $n = 0; $n < $ctrlnum; $n++){
			if(defined $windowtable{$chrom}{$key}{$sampname[$n]}) {
				$normfactor = $normfactor + $windowtable{$chrom}{$key}{$sampname[$n]}/$windowtable{"total"}{$sampname[$n]};
				if($windowtable{$chrom}{$key}{$sampname[$n]} < $mincoverage) {$countflag = 1;}
			}
			else{$countflag = 1;}
		}
		if($countflag == 1) {next;}
		$normfactor = (0.5 * $normfactor) / $ctrlnum;

		print OUT $chrom, "\t", $key*$windowsize, "\t", ($key+1)*$windowsize;
		for(my $n = 0; $n < @inputsam; $n++){
			my $cov = 0;
			if(defined $windowtable{$chrom}{$key}{$sampname[$n]}) {
				$cov = ($windowtable{$chrom}{$key}{$sampname[$n]} / $windowtable{"total"}{$sampname[$n]} ) / $normfactor;
				$cov = sprintf("%.2f", $cov);
			}
			print OUT "\t" , $cov;
		}
		print OUT "\n";
	}
}
close OUT;

__END__


#Calculate conversion rates
print "Calculating converian rates of all samples\n";
my %convertedtable;
my %allconvals;
foreach my $chrom (sort keys %windowtable){
	if($chrom eq "total") {next;}
	foreach my $key (sort { $a <=> $b } keys( %{$windowtable{$chrom}} ) ){
		my $countflag = 0;
		my $normfactor = 0;
		for(my $n = 0; $n < $ctrlnum; $n++){
			if(defined $windowtable{$chrom}{$key}{$sampname[$n]}) {
				$normfactor = $normfactor + $windowtable{$chrom}{$key}{$sampname[$n]}/$windowtable{"total"}{$sampname[$n]};
				if($windowtable{$chrom}{$key}{$sampname[$n]} < $mincoverage) {$countflag = 1;}
			}
			else{$countflag = 1;}
		}
		if($countflag == 1) {next;}
		$normfactor = (0.5 * $normfactor) / $ctrlnum;
		for(my $n = 0; $n < @inputsam; $n++){
			my $cov = 0;
			if(defined $windowtable{$chrom}{$key}{$sampname[$n]}) {
				$cov = ($windowtable{$chrom}{$key}{$sampname[$n]} / $windowtable{"total"}{$sampname[$n]} ) / $normfactor;
				$cov = sprintf("%.2f", $cov);
			}
			$convertedtable{$chrom}{$key}{$sampname[$n]} = $cov;
			push(@{$allconvals{$sampname[$n]}},$cov);
		}
	}
}
undef %windowtable;

# Get Medians
print "Calculating medians of all coverage\n";
my @median = ();
for(my $n = 0; $n < @inputsam; $n++){
	my $med = median( @{$allconvals{$sampname[$n]}});
	push(@median,$med);
}
undef %allconvals;


print "Printing table of best guess ploidy (rounded to nearest whole number)\n";
foreach my $chrom (sort keys %convertedtable){
	foreach my $key (sort { $a <=> $b } keys( %{$convertedtable{$chrom}} ) ){
		my %ploids = ();
		print OUT $chrom, "\t", $key*$windowsize, "\t", ($key+1)*$windowsize;
		for(my $n = 0; $n < @inputsam; $n++){
			my $cov = 0;
			if(defined $convertedtable{$chrom}{$key}{$sampname[$n]}) {
				$cov = (2 * $convertedtable{$chrom}{$key}{$sampname[$n]}) / $median[$n];
				$cov = floor($cov + .5);
				push(@{$ploids{$sampname[$n]}}, $cov);
			}
			if(@{$ploids{$sampname[$n]}} >= $smoother){
				my $smoothavg = 0;
				for(my $t = 0; $t < $smoother; $t++){
					$smoothavg = $smoothavg + $ploids{$sampname[$n]}[$t];
				}
				$smoothavg = floor(($smoothavg / $smoother) + .5);
				print OUT "\t" , $smoothavg;
				shift(@{$ploids{$sampname[$n]}});
			}
			else{
				print OUT "\t" , "2";
			}
		}
		print OUT "\n";
	}
}
close OUT;

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}





__END__
