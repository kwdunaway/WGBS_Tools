#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Last Update Date: 1-26-2016
# Version: 2.1
#
# Takes SAM output from BS_Seeker2 and creates percentage methylation BED files that
# can be uploaded to the UCSC genome browser or further analyzed through StochHMM.
#
# PCR duplicate filter: This script takes the first X reads as defined in parameter #7
#
# The positions in the resulting percent methylation (permeth) BED files are what you
# would get if you go to the following website. For example, if you go here: 
#    http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chrY:59032572,59032573
# if would return CG. However, when you look at the position on the genome browser, the
# color will only cover 1 base (the 2nd one).
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following parameters:
    1) Input sorted SAM file
    2) Output files prefix (folder and prefix)
    3) Bed track prefix
    4) UCSC genome version (ex: hg38)
    5) Methylation type (CG, CH)
    6) Strand (combined, positive, or negative)
    7) # max duplicate reads (ex: 1)
" unless @ARGV == 7;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outprefix = shift(@ARGV);
my $bedprefix = shift(@ARGV);
my $genome_version = shift(@ARGV);
my $meth_type = shift(@ARGV);

# Choose methylation type to search, validation checking
my @searchchars;
if ($meth_type eq "CG") { @searchchars = ("X","x");}
elsif($meth_type eq "CH") { @searchchars = ("Z","z", "Y","y");}
else { die "Methylation type $meth_type is not one of:  CG or CH\n\n";}

# Validation checking on strand type
my $strand_type = shift(@ARGV);
unless (($strand_type eq "combined") || ($strand_type eq "positive") || ($strand_type eq "negative")) {
    die "Strand type $strand_type is not one of: combined, positive, or negative\n\n";
}

my $max_dup_reads = shift(@ARGV); # Maximum duplicate reads allowed for analysis

# Global Variables
my $currentchrom = "Not Set Yet";
my %Methylation;
my %positions;
my $count = 0;

# Previous read info
my $prevstart = 0;
my $prevstrand = "+";
my $prevmethstring = "";
my $dupcount = 1;

# Columns for formatting SAM files
my $chrc = 2;
my $startc = 3;
my $methc = 14;
my $strandc = 11;



#############
# Main Loop #
#############

# Note
# In the hash, 1's signify a methylated hit and 0's signify an unmethylated hit

if ($meth_type eq "CG") # Run this process if CG
{
	print "Scanning for 'CG'\n";
	while(<IN>){
		my @line = split("\t", $_);
		
		if ($line[0] =~ /@/) {next;}
		# Take in chromosome
		my $chrom = $line[$chrc];
	
		# Takes out weird chromosomes like chr#_random and ect.
		if ($chrom =~ /_/) {next;}

		my $start = $line[$startc];
		my $methstring = substr $line[$methc], 5;
		my $strand = substr $line[$strandc], 5,1;

		# Skips read if only looking at positive or negative strand
		if($strand eq "-" && $strand_type eq "positive") {next;}
		if($strand eq "+" && $strand_type eq "negative") {next;}

		# If duplicate line, take information until $max_dup_reads limit hit
		if($prevstart == $start && $prevstrand eq $strand) {
			$dupcount++;
			if($dupcount > $max_dup_reads){next;}
		}
		else {$dupcount = 1;}

		# On next chromosome, so print current one if not set yet
		if($chrom ne $currentchrom){
			if($currentchrom ne "Not Set Yet"){
				Print_MethylationHash(\%Methylation, $outprefix, $currentchrom, $bedprefix);
			}
			# Reset variables
			%Methylation = ();
			$currentchrom = $chrom;
			$dupcount = 1;
			print "Starting " , $chrom , "\n";
		}
		
		# If found characters for CG (X/x), add to the hash
		if($prevmethstring =~ m/$searchchars[0]/ || $prevmethstring =~ m/$searchchars[1]/){
			Addto_MethylationHash(\%Methylation, $searchchars[0], $prevmethstring, $prevstart, $prevstrand);
		}
		# Increment
		$prevmethstring = $methstring;
		$prevstart = $start;
		$prevstrand = $strand;
	}
	# Finish last line/chromosome
	Addto_MethylationHash(\%Methylation, $searchchars[0], $prevmethstring, $prevstart, $prevstrand);
	Print_MethylationHash(\%Methylation, $outprefix, $currentchrom, $bedprefix);
}
elsif ($meth_type eq "CH") # Run this process instead if CH
{
	print "Scanning for 'CH'\n";
	while(<IN>){
		my @line = split("\t", $_);
		# Take in chromosome
		my $chrom = $line[$chrc];
	
		# Takes out weird chromosomes like chr#_random and ect.
		if ($chrom =~ /_/) {next;}

		my $start = $line[$startc];
		my $methstring = substr $line[$methc], 5;
		my $strand = substr $line[$strandc], 5,1;

		# Check for valid CH methylation characters
		next unless ($methstring =~ m/$searchchars[0]/ || $prevmethstring =~ m/$searchchars[1]/ || $methstring =~ m/$searchchars[2]/ || $methstring =~ m/$searchchars[3]/);

		# Skips read if only looking at positive or negative strand
		if($strand eq "-" && $strand_type eq "positive") {next;}
		if($strand eq "+" && $strand_type eq "negative") {next;}

		# If duplicate line, take longest read and then skip
		if($prevstart == $start && $prevstrand eq $strand) {
			if(length($methstring) > length($prevmethstring)){
				$prevmethstring = $methstring;
			}
			next;
		}

		# On next chromosome, so print current one if not set yet
		if($chrom ne $currentchrom){
			if($currentchrom ne "Not Set Yet"){
				Print_MethylationHash_CH(\%Methylation, $outprefix, $currentchrom, $bedprefix);
			}
			# Reset variables
			%Methylation = ();
			%positions = ();
			$currentchrom = $chrom;
			print "Starting " , $chrom , "\n";
		}

		# Get end position
		my $end = $line[5];
		$end =~ s/[^\d]//g;
		$end = $end + $start;
		my $positionstart = $start;
		# Flip start if negative strand to store strand type in hash
		if($strand eq "-") {$positionstart = $start * -1;}

		# If methylated, add to methylation hash
		if($prevmethstring =~ m/$searchchars[0]/ || $prevmethstring =~ m/$searchchars[2]/){
			Addto_MethylationHash_CH(\%Methylation, $prevmethstring, $prevstart, $prevstrand);
			$positions{$positionstart} = $end;
		}
		elsif($prevmethstring =~ m/$searchchars[1]/ || $prevmethstring =~ m/$searchchars[3]/){
			$positions{$positionstart} = $end;
		}

		# Increment
		$prevmethstring = $methstring;
		$prevstart = $start;
		$prevstrand = $strand;

		# Keep the hash small to not overload memory
		# Add to the hash every 1000 counts and empty positions hash
		$count++;
		if($count >= 1000)
		{
			$count = 0;
			foreach my $posstart (keys %positions)
			{
				if($positions{$posstart} < $start)
				{
					#print "Check: $posstart\t";
					#print "$positions{$posstart}\n";
					foreach my $posch (keys %Methylation)
					{
						if($posch < 0 && $posstart < 0)
						{
							if(abs($posch) >= abs($posstart)-1 && abs($posch) <= $positions{$posstart}-1)
							{
								$Methylation{$posch} = $Methylation{$posch} . "0";
								#print $posch,"\t",$Methylation{$posch},"\n";
							}
						}
						elsif($posch > 0 && $posstart > 0)
						{
							if(abs($posch) >= abs($posstart) && abs($posch) <= $positions{$posstart})
							{
								$Methylation{$posch} = $Methylation{$posch} . "0";
								#print $posch,"\t",$Methylation{$posch},"\n";
							}
						}
					}
					delete $positions{$posstart};
				}
			}
		}
	}
	# Last Time
	Addto_MethylationHash_CH(\%Methylation, $prevmethstring, $prevstart, $prevstrand);
	foreach my $posstart (keys %positions)
	{
		foreach my $posch (keys %Methylation)
		{
			#print $posch, "\n";
			if($posch < 0 && $posstart < 0)
			{
				if(abs($posch) >= abs($posstart)-1 && abs($posch) <= $positions{$posstart}-1)
				{
					$Methylation{$posch} = $Methylation{$posch} . "0";
					#print $posch,"\t",$Methylation{$posch},"\n";
				}
			}
			elsif($posch > 0 && $posstart > 0)
			{
				if(abs($posch) >= abs($posstart) && abs($posch) <= $positions{$posstart})
				{
					$Methylation{$posch} = $Methylation{$posch} . "0";
					#print $posch,"\t",$Methylation{$posch},"\n";
				}
			}
		}
		delete $positions{$posstart};
	}
	Print_MethylationHash_CH(\%Methylation, $outprefix, $currentchrom, $bedprefix);
}


###############
# Subroutines #
###############

# Add for CG
sub Addto_MethylationHash{
	my ($Methylation_ref, $charsearch, $methstring, $start, $strand) = @_;
	
	# Finds Methylated
	$charsearch =~ tr/a-z/A-Z/; # Search uppercases
	my $offset = 0;
	my $position = 0;
	# Search the string until it ends
	while ($position >= 0)
	{
  		$position = index($methstring, $charsearch, $offset);
		# If past string, done searching
  		if($position == -1) {last;}
		# If positive strand
 		if($strand eq "+"){
			my $startpos = $start + $position;
			# Exists already, append a 1 for methylated hit
			if(defined $Methylation_ref->{$startpos}) {$Methylation_ref->{$startpos} = $Methylation_ref->{$startpos} . "1";}
			# else, new find, make it 1 for a methylated hit
			else {$Methylation_ref->{$startpos} = "1";}
 		}
		# If negative strand
 		elsif($strand eq "-"){
			my $startpos = $start + length($methstring) - $position -2;
			if(defined $Methylation_ref->{$startpos}) {$Methylation_ref->{$startpos} = $Methylation_ref->{$startpos} . "1";}
			else {$Methylation_ref->{$startpos} = "1";}
 		}
		# else error, confused strand
 		else { die "Strand not + or - but $strand \n";}
		$offset=$position+1;
	}
	
	# Finds Unmethylated
	# Same algorithm as above but searches lowercase and adds 0's
	$charsearch =~ tr/A-Z/a-z/;
	$offset = 0;
	$position = 0;
	while ($position >= 0)
	{
  		$position = index($methstring, $charsearch, $offset);
  		if($position == -1) {last;}
 		if($strand eq "+"){
			my $startpos = $start + $position;
			if(defined $Methylation_ref->{$startpos}) {$Methylation_ref->{$startpos} = $Methylation_ref->{$startpos} . "0";}
			else {$Methylation_ref->{$startpos} = "0";}
 		}
 		elsif($strand eq "-"){
			my $startpos = $start + length($methstring) - $position -2;
			if(defined $Methylation_ref->{$startpos}) {$Methylation_ref->{$startpos} = $Methylation_ref->{$startpos} . "0";}
			else {$Methylation_ref->{$startpos} = "0";}
 		}
 		else { die "Strand not + or - but $strand \n";}
		$offset=$position+1;
	}
}

# Add for CH
# Similar algorithm as above for CG but only looks at methylated and at Z/Y
sub Addto_MethylationHash_CH{
	my ($Methylation_ref, $methstring, $start, $strand) = @_;

	my @search = ("Z", "Y");
	# Run for Z and Y
	foreach my $charsearch(@search)
	{
		my $offset = 0;
		my $position = 0;
		while ($position >= 0)
		{
  			$position = index($methstring, $charsearch, $offset);
  			if($position == -1) {last;}
 			if($strand eq "+"){
				my $startpos = $start + $position;
				if(defined $Methylation_ref->{$startpos}) {$Methylation_ref->{$startpos} = $Methylation_ref->{$startpos} . "1";}
				else {$Methylation_ref->{$startpos} = "1";}
 			}
 			elsif($strand eq "-"){
				my $startpos = -($start + length($methstring) - $position -2);
				if(defined $Methylation_ref->{$startpos}) {$Methylation_ref->{$startpos} = $Methylation_ref->{$startpos} . "1";}
				else {$Methylation_ref->{$startpos} = "1";}
 			}
 			else { die "Strand not + or - but $strand \n";}
			$offset=$position+1;
		}
	}
}

# Print the hash to the output file for CG methylation
sub Print_MethylationHash{
	my ($Methylation_ref, $outprefix, $currentchrom, $bedprefix) = @_;
	
	my $outfile = $outprefix . "_" . $currentchrom . ".bed";
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	# Print header
	print OUT "track name=" , $bedprefix, $currentchrom, " description=" , $bedprefix, "_", $currentchrom, " useScore=0 itemRgb=On db=" , $genome_version , "\n";
	# For each position
	foreach my $posstart (sort { $a <=> $b }  keys %{$Methylation_ref}) {
		# Example print format
		#chr10   51332   51333   0.50-2  0       +       0       0       27,74,210
		my $posend = $posstart + 1;
		# Calculate the percentage methylation
		my $methperc = 0;
		my @methraw = split("",$Methylation_ref->{$posstart});
		my $methnum = @methraw;
		while(@methraw){
			$methperc += $methraw[0];
			shift(@methraw);
		}
		$methperc = $methperc / $methnum;
		$methperc = sprintf("%.2f", $methperc);
		# Choose the color, low meth > blue, med > green, high > red
		my $color = "0,0,0"; #black
		if ($methperc > 0 && $methperc <= .6) {$color = "27,74,210";} #blue
		elsif ($methperc > .6 && $methperc <= .8) {$color = "27,210,57";} #green
		elsif ($methperc > .8) {$color = "210,27,27";} #red
		# Finally print
		print OUT $currentchrom , "\t" , $posstart , "\t" , $posend , "\t" , $methperc , "-", $methnum , "\t" , "0\t+\t0\t0\t" , $color , "\n";
	}
	close OUT;
}

# Print the hash to the output file for CH methylation
# Similar algorithm to the above 
sub Print_MethylationHash_CH{
	my ($Methylation_ref, $outprefix, $currentchrom, $bedprefix) = @_;
	
	my $outfile = $outprefix . "_" . $currentchrom . ".bed";
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	print OUT "track name=" , $bedprefix, $currentchrom, " description=" , $bedprefix, "_", $currentchrom, " useScore=0 itemRgb=On db=" , $genome_version , "\n";
	foreach my $posstart (sort { abs($a) <=> abs($b) }  keys %{$Methylation_ref}) {
		#Example print format
		#chr10   51332   51333   0.50-2  0       +       0       0       27,74,210
		my $str = "+";
		# Find the start (depends on strand)
		my $trueposstart = $posstart;
		if($posstart < 0)
		{
			$trueposstart = ($posstart * -1);
			$str = "-";
		}
		my $posend = $trueposstart + 1;
		# Calculate percentage methylation
		my $methperc = 0;
		my @methraw = split("",$Methylation_ref->{$posstart});
		my $methnum = @methraw;
		foreach my $falsezero (@methraw)
		{
			if($falsezero == 1) {$methnum--;}
		}
		while(@methraw){
			$methperc += $methraw[0];
			shift(@methraw);
		}
		$methperc = $methperc / $methnum;
		$methperc = sprintf("%.2f", $methperc);
		# Choose the color, negative strand is yellow, positive is magenta
		my $color = "0,0,0"; #black
		if ($methperc > 0 && $str eq "-") {$color = "0,51,0";} #yellow
		elsif ($methperc > 0 && $str eq "+") {$color = "255,0,255";} #magenta
		# Finally print if methylation percentage is not 0
		if ($methperc != 0) 
		{
			print OUT $currentchrom , "\t" , $trueposstart , "\t" , $posend , "\t" , $methperc , "-", $methnum , "\t" , "0\t", $str, "\t0\t0\t" , $color , "\n";
		}
	}
	close OUT;
}
