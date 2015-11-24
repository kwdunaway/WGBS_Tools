#!/usr/bin/perl 
use strict; use warnings; 
use POSIX 'floor';
use Time::HiRes;
##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu and rgchu@ucdavis.edu
# Date: 2-24-2015
#
# *INSTRUCTIONS*
# Run this script (window_percentage_methylation.pl) first
# And then the bash script (makeMyWindows.bash)
#
# This script takes windows (user defined parameters) and outputs avg meth across windows
# based on a read centric method. The script also outputs a count of CpG assays.
#
# Note: This script produces multiple temporary files that are deleted
# automatically after use. These files are the temporary bash file
# (named makeMyWindows.bash) and as many temporary chromosome files as
# exist in the input ((output name)+"chr"+(chromosome number))
#
# For example: If 2 CpGs were assayed like this:
#   0.5-2
#   1-8
# The output is .75 for methylation and 2 for coverage.
# This script would output .9 for methylation and 10 for coverage.
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################
my $start = Time::HiRes::gettimeofday();
die "$0 usage, needs the following parameters: 
    1) Output table file name
    2) Window size (ex: 20000)
    3) Min # of CpGs per window (ex: 20)
    4) Min # of reads per CpG counted (ex: 1)
    5) Min # of samples that must meet requirements 3 and 4
    6,8+) Permeth prefix (ex: JLKD001/PerMeth_JLKD001/PerMeth_JLKD001_chr)
    7,9+) Name of experiments in output file (ex: JLKD001)

    Note: must run this script while in the same folder as the script
" unless @ARGV > 6;

my $outputname = shift(@ARGV);
my $windowsize = shift(@ARGV);
my $mincpg = shift(@ARGV);
my $mincoverage = shift(@ARGV);
my $minfiles = shift(@ARGV);
my @Permethfiles; 
my @Permethnames; 
while(@ARGV){
	push(@Permethfiles, shift(@ARGV));
	push(@Permethnames, shift(@ARGV));
}

# Label output table header
open(OUT, ">$outputname") or die "Error: cannot open $outputname OUT outfile";
my $countoutputname = $outputname . ".count";
open(COUNT, ">$countoutputname") or die "Error: cannot open $countoutputname COUNT outfile";

my @samplenames;
print OUT "chr" , "\t" , "start" , "\t" , "end";
print COUNT "chr" , "\t" , "start" , "\t" , "end";
for (my $n = 0; $n < @Permethnames; $n++){
	print OUT "\t" , $Permethnames[$n];
	print COUNT "\t" , $Permethnames[$n];
}
print OUT "\n";
print COUNT "\n";

close OUT;
close COUNT;

# Construct temporary bash to run scripts in parallel and speed up run time
my $outputbash = "./makeMyWindows.bash";
open(BASH, ">$outputbash") or die "$0: Error: cannot open $outputbash BASH outfile";

# Scan directory for number of chromosomes
my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
my $filedir = $Permethfiles[0];
$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
my @files = glob( $filedir . '*' ); # Gets list of all files in directory
@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

foreach my $file (@files) {
	$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
	push (@Chr, $1); # Add to list of chromosome numbers
}
@Chr = sort @Chr; # Sort list of chromosome numbers

# Run for each chromosome (& runs in background simultaneously)
my $processcurrent = 0; # Control the number of processes run simultaneously
my $processlimit = 2; # Change this number to change the limit on processes
foreach my $chr (@Chr){
	print BASH "perl ./window_permeth_onechr.pl ", $outputname, " ", $chr, " ", $windowsize, " ", $mincpg, " ", $mincoverage, " ", $minfiles;
	for(my $i = 0; $i < @Permethfiles; $i++){
		print BASH " ", $Permethfiles[$i], " ", $Permethnames[$i];
	}
	print BASH " &\n";
	$processcurrent++;
	if($processcurrent == $processlimit){
		print BASH "wait\n";
		$processcurrent = 0;
	}
}

# Wait until all processes are finished
print BASH "\nwait\n";
print BASH "echo Finished\n";

# Append output files to the aggregate output table
print BASH "cat $outputname ";
foreach my $chr (@Chr){
	print BASH $outputname, "chr", $chr, " ";
}
print BASH "> ", $outputname, ".temp\n";
print BASH "mv $outputname.temp $outputname\n";

print BASH "cat $countoutputname ";
foreach my $chr (@Chr){
	print BASH $outputname, "chr", $chr, ".count ";
}
print BASH "> ", $countoutputname, ".temp\n";
print BASH "mv $countoutputname.temp $countoutputname\n";

# Remove single files
foreach my $chr (@Chr){
	print BASH "rm ", $outputname, "chr", $chr, "\n";
	print BASH "rm ", $outputname, "chr", $chr, ".count", "\n";
}
print BASH "echo Compiled\n";

# Self-delete temporary bash file
print BASH "rm -- \"\$0\"\n";

close BASH;

# Give permissions to run bash file
system('chmod 777 ./makeMyWindows.bash');

# Run bash file (which deletes upon completion)
system('./makeMyWindows.bash');

# Get run time
my $end = Time::HiRes::gettimeofday();
print "TIMERPL FINISH\n";
printf("%.2f\n", $end - $start);
