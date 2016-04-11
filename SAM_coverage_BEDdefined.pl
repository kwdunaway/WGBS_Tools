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
    2) Bed file to define areas
    3+) Input SAM File(s)
" unless @ARGV > 2;

my $outputfile = shift(@ARGV); # Name of output file
open(OUT, ">$outputfile") or die "I/O Error: cannot open $outputfile output file";
my $inputfile = shift(@ARGV); # Name of output file
open(IN, "<$inputfile") or die "I/O Error: cannot open $inputfile output file";
my @inputsam = @ARGV;

# Global variables
my %windowtable;
# Format windowtable{chr}{start}{sample} = count;
my $filecount = scalar @ARGV;
my $sizeinputsam = $#inputsam+1;
my @sampname = ();

# Header
print "Initializing.\n";
my $header = "chr\tstart\tend\tname";
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
print OUT "Total\tTotal\tTotal\tTotal";



#####################################
#  Input bed file for defined areas #
#####################################

my $namecount = 1;
my $firstline = <IN>; 
my @line = split("\t",$firstline);
if($line[0] =~ /^chr/){
	print "Header not detected in first line, processing it.\n";
	my $name = $namecount;
	if(scalar @line > 3){$name = $line[3]}
	else{$namecount++;}
	$windowtable{$line[0]}{$line[1]}{$line[2]}{"name"} = $name;	
}
else {print "Header detected in first line, skipping it.\n";}

while(<IN>){
	@line = split("\t",$_);
	my $name = $namecount;
	if(scalar @line > 3){$name = $line[3]}
	else{$namecount++;}
	$windowtable{$line[0]}{$line[1]}{$line[2]}{"name"} = $name;
#	print "before\n";
	for(my $n = 0; $n < @sampname; $n++){
#		print "for loop\n";
		$windowtable{$line[0]}{$line[1]}{$line[2]}{$sampname[$n]} = 0;
	}
}



###############
#  Processing #
###############

for(my $n = 0; $n < @inputsam; $n++){
	$inputfile = $inputsam[$n];
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
		
		print "";
		foreach my $keystart (sort { $a <=> $b } keys( %{$windowtable{$chrom}} ) ){
			if($start > $keystart){
				foreach my $keyend (sort { $a <=> $b } keys( %{$windowtable{$chrom}{$keystart}} ) ){
					if($start < $keyend){
						$windowtable{$chrom}{$keystart}{$keyend}{$sampname[$n]}++;
					}
				}
			}
		}

	}
	close SAM;
	print OUT "\t" , $windowtable{"total"}{$sampname[$n]};
}
print OUT "\n";



#######################
#  Printing OUT table #
#######################

print "Printing OUT table\n";
foreach my $chrom (sort keys %windowtable){
	if($chrom eq "total") {next;}
	foreach my $keystart (sort { $a <=> $b } keys( %{$windowtable{$chrom}} ) ){
		foreach my $keyend (sort { $a <=> $b } keys( %{$windowtable{$chrom}{$keystart}} ) ){
			print OUT $chrom , "\t" , $keystart , "\t" , $keyend , "\t" , $windowtable{$line[0]}{$line[1]}{$line[2]}{"name"};
			for(my $n = 0; $n < @sampname; $n++){
				print OUT "\t" , $windowtable{$line[0]}{$line[1]}{$line[2]}{$sampname[$n]};
			}
			print OUT "\n";			
		}
	}
}

__END__


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
