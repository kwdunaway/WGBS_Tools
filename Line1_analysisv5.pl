#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Version: 4.0
# Last Updated: 1-6-2015
#
# This script looks through all raw fastq sequencing reads and finds the reads that
# have the Line1 pattern.  Then, it can put those reads in a separate file.  Also, the
# script quantifies methylation of these sequences across four potential 
# methylation sites.
#
# To call:
#   perl /data/scratch/programs/perl_script/Line1_analysisv2.pl YES Stats_L1v2.txt Sample_JLD018/JLD018_filtered.fq_LINE1reads.fq Sample_JLDS019/JLD019_filtered.fq_LINE1reads.fq Sample_JLKD006/JLKD_006_filtered.fq_LINE1reads.fq Sample_JLD017/JLD017_filtered.fq_LINE1reads.fq Sample_JLKD008/JLKD008_filtered.fq_LINE1reads.fq Sample_JLKD007/JLKD007_filtered.fq_LINE1reads.fq Sample_JLKD002/JLKD002_filtered.fq_LINE1reads.fq Sample_JLKD003/JLKD003_filtered.fq_LINE1reads.fq Sample_JLKD005/JLKD_005_filtered.fq_LINE1reads.fq Sample_JLKD001/JLKD001_filtered.fq_LINE1reads.fq Sample_JLKD004/JLKD_004_filtered.fq_LINE1reads.fq
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Filter (If NO, nothing. If YES, creates a Line-1 fastq file with line1.fq added to the end)
    2) Results Table Outfile 
    3-?) Input fastq file(s) (only uses first if Filtering)
" unless @ARGV > 2;

my $filtertag = shift(@ARGV);	
$filtertag =~ tr/[a-z]/[A-Z]/;
my $results_outfile = shift(@ARGV);	
open(RESULTS, ">$results_outfile") or die "cannot open $results_outfile outfile";
my @infiles = @ARGV;

# Hard coded for now, but should be able to be set as parameter
my $maxcopies = 1;

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

print RESULTS "Sample\tcount\tMeth1\tMeth2\tMeth3\tMeth4\n";
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
	if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
	else{open(IN, "<$infile") or die "cannot open $infile infile";}
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
		$printstring = $printstring . "\t" . $SNPsequences{$searchstring}{"count"};
		if($SNPsequences{$searchstring}{"count"} == 0){next;}
		for(my $SNPloc = 1; $SNPloc < 5; $SNPloc++){
			$printstring = $printstring . "\t" . sprintf("%.3f", $SNPsequences{$searchstring}{$SNPloc}{"C"}/$SNPsequences{$searchstring}{"count"});
		}
		$printstring = $printstring . "\n";
	}
	return($printstring);
}

__END__
sub BS_Convert_All{
    my $dna = shift;
    
    #Convert everything to uppercase
    $dna =~ tr/[a-z]/[A-Z]/;
	#Get reverse complement
    my $rcdna = reverse_complement($dna);
    
    #Hash with multiple information about sequence to assay in it
    my %sequences;

    #Get bisulfite converted sequence, both strand conversion
    $sequences{"start"}{"seq"} = $dna;
    $sequences{"revcomp"}{"seq"} = $rcdna;
    $sequences{"bsfs"}{"seq"} = bs_convert_same($dna);
    $sequences{"bsfo"}{"seq"} = bs_convert_opposite($dna);
    $sequences{"bsrs"}{"seq"} = bs_convert_same($rcdna);
    $sequences{"bsro"}{"seq"} = bs_convert_opposite($rcdna);

	#Get search string for CpG's
    $sequences{"bsfs"}{"search"} = annotate_cpgs($sequences{"bsfs"}{"seq"});
    $sequences{"bsfo"}{"search"} = annotate_cpgs($sequences{"bsfo"}{"seq"});
    $sequences{"bsrs"}{"search"} = annotate_cpgs($sequences{"bsrs"}{"seq"});
    $sequences{"bsro"}{"search"} = annotate_cpgs($sequences{"bsro"}{"seq"});
    
    #Find SNP sites
    @{$sequences{"bsfs"}{"SNPs"}} = find_SNPsites($sequences{"bsfs"}{"seq"}, 'Y');
    @{$sequences{"bsfo"}{"SNPs"}} = find_SNPsites($sequences{"bsfo"}{"seq"}, 'R');
    @{$sequences{"bsrs"}{"SNPs"}} = find_SNPsites($sequences{"bsrs"}{"seq"}, 'Y');
    @{$sequences{"bsro"}{"SNPs"}} = find_SNPsites($sequences{"bsro"}{"seq"}, 'R');

#    my %sequences = ($dna_bsfs,$dna_bsfo,$dna_bsrs,$dna_bsro);
	return(%sequences);
}

sub reverse_complement {
    my $dna = shift;

	# reverse the DNA sequence
    my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

    return $revcomp;
}

sub bs_convert_same {
    my $bsdna = shift;

    #Put SNP Y (C or T) for every CpG (because the C's might be methylated
    $bsdna =~ s/CG/YG/g;

    #Convert all non CpGs C's to T's (because they won't be methylated)
    $bsdna =~ s/C/T/g;

    return $bsdna;
}

sub bs_convert_opposite {
    my $bsdna = shift;

    #Put SNP R (G or A) for every CpG (because the C's might be methylated on opposite strand
    $bsdna =~ s/CG/CR/g;

    #Convert all non CpGs G's to A's (because the C's won't be methylated)
    $bsdna =~ s/G/A/g;

    return $bsdna;
}

#Returns SNP locations as an array of numbers
sub find_SNPsites{
	# String containing sequence to be looked through for a SNP
	my $string = shift;
	# Character to be looked through the string for, most likely Y for R.
	my $char = shift;
	# Resulting locations to be returned
	my @locations;

	my $offset = 0;
	my $result = index($string, $char, $offset);
	while ($result != -1) {
		push(@locations, $result);
		$offset = $result + 1;
    	$result = index($string, $char, $offset);
  }
	
	return(@locations);
}

sub annotate_cpgs{
    my $searchstring = shift;
    $searchstring =~ s/Y/[CT]/g;
    $searchstring =~ s/R/[GA]/g;
    $searchstring = "(" . $searchstring . ")";
    return($searchstring);
}

sub Print_BS_Sequences {
	my ($seq_ref) = @_;
#    my $params = shift;
    my %sequences = %$seq_ref;
    
	my $print_statement =
	"Original:             \t" . $sequences{"start"}{"seq"} . "\n" .
	"BS For Same:          \t" . $sequences{"bsfs"}{"seq"} . "\n" .
	"Searchstring:         \t" . $sequences{"bsfs"}{"search"} . "\n\n" .
	
	"Original:             \t" . $sequences{"start"}{"seq"} . "\n" .
	"BS For Opp:           \t" . $sequences{"bsfo"}{"seq"} . "\n" .
	"Searchstring:         \t" . $sequences{"bsfo"}{"search"} . "\n\n" .

	"Reverse:              \t" . $sequences{"revcomp"}{"seq"} . "\n" .
	"BS Rev Same:          \t" . $sequences{"bsrs"}{"seq"} . "\n" .
	"Searchstring:         \t" . $sequences{"bsrs"}{"search"} . "\n\n" .
	
	"Reverse:              \t" . $sequences{"revcomp"}{"seq"} . "\n" .
	"BS Rev Opp:           \t" . $sequences{"bsro"}{"seq"} . "\n" . 
	"Searchstring:         \t" . $sequences{"bsro"}{"search"} . "\n\n";
	
	return($print_statement);
}

sub GetAllSNP {
    my $infile = shift;
    my $maxcopies = shift;
	if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
	else{open(IN, "<$infile") or die "cannot open $infile infile";}
	my ($seq_ref) = @_;
    my %sequences = %$seq_ref;
    my %SNPhash;
    
    my $counter = 0;

	while (<IN>) {
		$counter++;
		if($counter % 100000 == 0) {print $counter , "\n";}
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;	    

		foreach my $key (keys %sequences) {
	    	if($seq =~ /$sequences{$key}{"search"}/){
				my $match = $1;
				my $searchterm = $sequences{$key}{"search"};
				if (defined $SNPhash{$seq}){
					$SNPhash{$key}{$seq}{"copy_number"}++;
				}
				else{
					$SNPhash{$key}{$seq}{"copy_number"} = 1;
					@{$SNPhash{$key}{$seq}{"SNPs"}} = FindSeqSNPs($searchterm, @{$sequences{$key}{"SNPs"}});
				}
			}
		}
	}
	close IN;
	return(%SNPhash);
}

sub GetAllMatchReads {
    my $infile = shift;
	if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
	else{open(IN, "<$infile") or die "cannot open $infile infile";}
	my @sequences_array = @_;
    my @matched_reads;
    my @counts;
    for (my $r = 0; $r < @sequences_array; $r++){
		$matched_reads[$r] = "";
		$counts[$r] = 0;
    }
    my $counter = 0;
	while (<IN>) {
		$counter++;
		if($counter % 100000 == 0) {print $counter , "\n";}
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;
		for(my $n = 0; $n < @sequences_array; $n++){
	    	if($seq =~ /$sequences_array[$n]/){
	    		$matched_reads[$n] = $matched_reads[$n] . $ID . $seq . "\n" . $third . $quality;
	    		$counts[$n]++;
			}
		}
	}
	close IN;
	
	for(my $n = 0; $n < @sequences_array; $n++){ print "$sequences_array[$n]\t\t$counts[$n]\n"; }
	
	return(@matched_reads);
}

sub CaptureSNPs {
    my $infile = shift;
	if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
	else{open(IN, "<$infile") or die "cannot open $infile infile";}
	my @sequences_array = @_;
    my %SNPsequences;
#    my @counts;
#    for (my $r = 0; $r < @sequences_array; $r++){
#		$matched_reads[$r] = "";
#		$counts[$r] = 0;
#    }
#    my $counter = 0;
	while (<IN>) {
#		$counter++;
#		if($counter % 100000 == 0) {print $counter , "\n";}
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;
		for(my $n = 0; $n < @sequences_array; $n++){
	    	if($seq =~ /$sequences_array[$n]/){
	    		if(defined $SNPsequences{$sequences_array[$n]}{$1}){ $SNPsequences{$sequences_array[$n]}{$1}++; }
	    		else{ $SNPsequences{$sequences_array[$n]}{$1} = 1; }
#	    		$matched_reads[$n] = $matched_reads[$n] . $ID . $seq . "\n" . $third . $quality;
#	    		$counts[$n]++;
			}
		}
	}
	close IN;
	
#	for(my $n = 0; $n < @sequences_array; $n++){ print "$sequences_array[$n]\t\t$counts[$n]\n"; }
	
	return(%SNPsequences);
}

# returns a hash of characters at specific location in a string
sub FindSeqSNPs {
	my $string = shift;
	my @SNPsites = @_;
	my @Results;
	while(@SNPsites){
		my $char = shift(@SNPsites);
		push (@Results ,substr($string, $char, 1));
	}
	return(@Results);
}
	
sub ConsolidateBSSNPs{
	my ($seq_ref) = @_;
    my %SNPhash = %$seq_ref;
    #my @names = ("bsfs","bsfo","bsrs","bsro");
    my %running_results = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsfs"}});
    my %bsfs = %running_results;
    
#    print "bsfs           \tT=" , $bsfs{"0"}{"T"}, " C=" , $bsfs{"0"}{"C"} , "\tT=", $bsfs{"1"}{"T"}, " C=" , $bsfs{"1"}{"C"} , "\tT=", $bsfs{"2"}{"T"}, " C=" , $bsfs{"2"}{"C"} , "\tT=", $bsfs{"3"}{"T"}, " C=" , $bsfs{"3"}{"C"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";
    
    #BS conversion in forward orientation on opposite strand
    my %bsfo = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsfo"}});
	my $n = -1;
    while(1){
    	$n++;
    	if(defined $running_results{$n}){
	   		if(defined $bsfo{$n}{"C"} or $bsfo{$n}{"T"}){
				foreach my $key (keys %{$bsfo{$n}}) {
				    if(defined $running_results{$n}{$key}){
		    			$running_results{$n}{$key} = $running_results{$n}{$key} + $bsfo{$n}{$key};
	    			}
	   				else {$running_results{$n}{$key} = $bsfo{$n}{$key};}
				}
			}
    		else{
		    	if(defined $running_results{$n}{"C"} && defined $bsfo{$n}{"G"}){
	    			$running_results{$n}{"C"} = $running_results{$n}{"C"} + $bsfo{$n}{"G"};
	    		}
	    		else {$running_results{$n}{"C"} = $bsfo{$n}{"G"};}
		    	if(defined $running_results{$n}{"T"} && defined $bsfo{$n}{"A"}){
	   	 			$running_results{$n}{"T"} = $running_results{$n}{"T"} + $bsfo{$n}{"A"};
	    		}
	    		else {$running_results{$n}{"T"} = $bsfo{$n}{"A"};}
	    	}
	    }
    	else{last;} 
    }
    my $SNPnum = $n;
#    print "bsfo           \tT=" , $bsfo{"0"}{"A"}, " C=" , $bsfo{"0"}{"G"} , "\tT=", $bsfo{"1"}{"A"}, " C=" , $bsfo{"1"}{"G"} , "\tT=", $bsfo{"2"}{"A"}, " C=" , $bsfo{"2"}{"G"} , "\tT=", $bsfo{"3"}{"A"}, " C=" , $bsfo{"3"}{"G"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";

    #BS conversion in reverse orientation on same strand 
    my %bsrs = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsrs"}});
	$n = -1;
	my $t = $SNPnum;
    while(1){
    	$n++;
    	$t--;
    	if(defined $running_results{$t}){
    		if(defined $bsrs{$n}{"A"} or $bsrs{$n}{"G"}){
				foreach my $key (keys %{$bsrs{$n}}) {
			    	if(defined $running_results{$t}{$key}){
		    			$running_results{$t}{$key} = $running_results{$t}{$key} + $bsrs{$n}{$key};
	    			}
	    			else {$running_results{$t}{$key} = $bsrs{$n}{$key};}
				}
   			}
    		else{
		    	if(defined $running_results{$t}{"C"} && defined $bsrs{$n}{"C"}){
	    			$running_results{$t}{"C"} = $running_results{$t}{"C"} + $bsrs{$n}{"C"};
	    		}
	    		else {$running_results{$t}{"C"} = $bsrs{$n}{"C"};}
		    	if(defined $running_results{$t}{"T"} && defined $bsrs{$n}{"T"}){
		    		print $running_results{$t}{"T"} , "\n";
		    		print $bsrs{$n}{"T"} , "\n";
	   	 			$running_results{$t}{"T"} = $running_results{$t}{"T"} + $bsrs{$n}{"T"};
	    		}
	    		else {$running_results{$t}{"T"} = $bsrs{$n}{"T"};}
	    	}
	    }
    	else{last;} 
    }
#    print "bsrs           \tT=" , $bsrs{"3"}{"T"}, " C=" , $bsrs{"3"}{"C"} , "\tT=", $bsrs{"2"}{"T"}, " C=" , $bsrs{"2"}{"C"} , "\tT=", $bsrs{"1"}{"T"}, " C=" , $bsrs{"1"}{"C"} , "\tT=", $bsrs{"0"}{"T"}, " C=" , $bsrs{"0"}{"C"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";
    
    #BS conversion in reverse orientation on opposite strand 
    my %bsro = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsro"}});
	$n = -1;
	$t = $SNPnum;
    while(1){
    	$n++;
    	$t--;
    	if(defined $running_results{$t}){
    		if(defined $bsro{$n}{"C"} or $bsro{$n}{"T"}){
				foreach my $key (keys %{$bsro{$n}}) {
			    	if(defined $running_results{$t}{$key}){
		    			$running_results{$t}{$key} = $running_results{$t}{$key} + $bsro{$n}{$key};
	    			}
	    			else {$running_results{$t}{$key} = $bsro{$n}{$key};}
				}
   		}
    		else{
		    	if(defined $running_results{$t}{"C"} && defined $bsro{$n}{"G"}){
	    			$running_results{$t}{"C"} = $running_results{$t}{"C"} + $bsro{$n}{"G"};
	    		}
	    		else {$running_results{$t}{"C"} = $bsro{$n}{"G"};}
		    	if(defined $running_results{$t}{"T"} && defined $bsro{$n}{"A"}){
	   	 			$running_results{$t}{"T"} = $running_results{$t}{"T"} + $bsro{$n}{"A"};
	    		}
	    		else {$running_results{$t}{"T"} = $bsro{$n}{"A"};}
	    	}
	    }
   	else{last;} 
    }
#    print "bsfo           \tT=" , $bsrs{"3"}{"A"}, " C=" , $bsrs{"3"}{"G"} , "\tT=", $bsrs{"2"}{"A"}, " C=" , $bsrs{"2"}{"G"} , "\tT=", $bsrs{"1"}{"A"}, " C=" , $bsrs{"1"}{"G"} , "\tT=", $bsrs{"0"}{"A"}, " C=" , $bsrs{"0"}{"G"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";

    my %results;
    %{$results{"total"}} = %running_results;
    %{$results{"bsfs"}} = %bsfs;
    %{$results{"bsfo"}} = %bsfo;
    %{$results{"bsrs"}} = %bsrs;
    %{$results{"bsro"}} = %bsro;
	return(%results);
}

# Returns a hash with the format SNPresults{numberSNPposition}{SNPletter}=count
sub ConsolidateSNPs{
	#   Example of format:
	#   SNPresults{0}{C} = 12
	#   SNPresults{0}{T} = 108
	# This would be 10% C and 90% T

	my $maxcopies = shift;
	my ($seq_ref) = @_;
    my %Readhash = %$seq_ref;
    my %SNPresults;

	foreach my $key (keys %Readhash) {

#		print $Readhash{$key}{"copy_number"} , "\n";

#		my $copiesused = $Readhash{$key}{"copy_number"};
#		if($copiesused > $maxcopies){
			my $copiesused = $maxcopies;
#		}
		my @tmparray = @{$Readhash{$key}{"SNPs"}};
		for(my $n = 0; $n < @tmparray; $n++){
			if(defined $SNPresults{$n}{$tmparray[$n]}){
				$SNPresults{$n}{$tmparray[$n]} = $SNPresults{$n}{$tmparray[$n]} + $copiesused;		
			}
			else{
				$SNPresults{$n}{$tmparray[$n]} = $copiesused;					
			}
		}
	}
	return(%SNPresults);
}


sub PrintSNPs{
	my ($seq_ref) = @_;
    my %results = %$seq_ref;
    
    my $printline = "";

#    my %printresults;    
#	foreach my $file (keys %results) {
#		foreach my $type (keys %{results{$file}}) {
#			%{$printresults{$type}{$file}} = %{$results{$file}{$type}}
#		}
#   }

	foreach my $file (keys %results) {
		foreach my $type (keys %{$results{$file}}) {
			$printline  = $printline . $file . "\t" . $type . "\t";
			foreach my $SNP (keys %{$results{$file}{$type}}) {
			   $printline = $printline . $SNP . "=" . $results{$file}{$type}{$SNP} . " ";
			 }
			 $printline = $printline . "\n";
		}
		$printline = $printline . "\n";
   	}    
	return($printline);
}


__END__

#################
# In Files Loop #
#################

while(@Infiles){
	my $infile = shift(@Infiles);
	open(IN, "<$infile") or die "cannot open $infile infile";
	print "\nStarting file: $infile\n";
#	print "Number or Reads matched so far:\t";

#	my $output_filename = $infile . $out_suffix;
#	open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";

###############################
# Initialization of Variables #
###############################

	my $meth1 = 0;
	my $meth2 = 0;
	my $meth3 = 0;
	my $meth4 = 0;
	my $unmeth1 = 0;
	my $unmeth2 = 0;
	my $unmeth3 = 0;
	my $unmeth4 = 0;
	my $SNPG = 0;
	my $SNPT = 0;
	my $LINE1 = 0;
	
	my %Reads;
	my %FlankBefore;
#	my %FlankAfter;
#	my %FlankBoth;

#############################
# Main Search Function Loop #
#############################

	while (<IN>)
	{
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;
	    if($seq =~ /(TT([CT])GTGGTG([CT])GT([CT])GTTTTTTAA([GT])T([CT])GGTT)/){
			if (defined $Reads{$seq}){
				$Reads{$seq}{"copy_number"}++;
			}
			else{
				# Initialize hash variables
				$Reads{$seq}{"copy_number"} = 1;
				$Reads{$seq}{"methylation"} = 0;
				$Reads{$seq}{"site1"}=0;
				$Reads{$seq}{"site2"}=0;
				$Reads{$seq}{"site3"}=0;
				$Reads{$seq}{"site4"}=0;
				
				# Get site methylation information data
				#Variable that counts and bins read methylation, 
				my $read_meth = 0;
				if($2 eq "C") {$meth1++; $read_meth++; $Reads{$seq}{"site1"}=1;} elsif($2 eq "T") {$unmeth1++;} else{die "Site1: $2";}
				if($3 eq "C") {$meth2++; $read_meth++; $Reads{$seq}{"site2"}=1;} elsif($3 eq "T") {$unmeth2++;} else{die "Site2: $3";}
				if($4 eq "C") {$meth3++; $read_meth++; $Reads{$seq}{"site3"}=1;} elsif($4 eq "T") {$unmeth3++;} else{die "Site3: $4";}
				if($6 eq "C") {$meth4++; $read_meth++; $Reads{$seq}{"site4"}=1;} elsif($6 eq "T") {$unmeth4++;} else{die "Site4: $6";}
				if($5 eq "G") {$SNPG++;} elsif($5 eq "T") {$SNPT++;} else{die "A/T: $5";}
				$Reads{$seq}{"methylation"} = $read_meth;

				# Get flanking region information
				my @Flank = split(/TT[CT]GTGGTG[CT]GT[CT]GTTTTTTAA[GT]T[CT]GGTT/, $seq, 2);
#				print $seq ,"\n"; 				
#				print $Flank[0], "                             " , $Flank[1], "\n";
				my $FlankUp = substr($Flank[0], -30);
#				my $FlankDown = substr($Flank[1], 0, 20);
#				print "Up:\t$FlankUp\nDown:\t$FlankDown\n\n";
				if(length($FlankUp) == 30){
					$FlankUp = substr($FlankUp, 0, 20);
					if (defined $FlankBefore{$read_meth}{$FlankUp}){$FlankBefore{$read_meth}{$FlankUp}++;} else {$FlankBefore{$read_meth}{$FlankUp}=1;}
				}
#				if(length($FlankDown) == 20){
#					if (defined $FlankAfter{$read_meth}{$FlankDown}){$FlankAfter{$read_meth}{$FlankDown}++;} else {$FlankAfter{$read_meth}{$FlankDown}=1;}
#				}
#				if (length($FlankUp) + length($FlankDown) == 40){
#					my $flankb = $FlankUp . "  " . $FlankDown;
#					if (defined $FlankBoth{$read_meth}{$flankb}) {$FlankBoth{$read_meth}{$flankb}++;} else {$FlankBoth{$read_meth}{$flankb}=1;}				
#				}
				
				$LINE1++;
#	    		if($LINE1 % 20 == 0) {print $LINE1 , "\n";}
			}
#			print OUT $ID , $seq , "\n", $third , $quality;
    	}
		else { die "$seq does not match up";}
	}
	
	
	for(my $read_meth = 0; $read_meth < 5; $read_meth++){
#		last;
#		next;
#		delete $FlankBefore{$read_meth}{"GGATTTTTTGAGTTAGGTGT"};

		# Upstream flanking region
		my $first = 0;	
		foreach my $key (sort { $FlankBefore{$read_meth}{$b} <=> $FlankBefore{$read_meth}{$a} } keys %{ $FlankBefore{$read_meth} }) {
			if($first eq "0") { 
				$first = $key; 
				print "Methylation:\t" , $read_meth, "\n";
				print $key , "\t" , $FlankBefore{$read_meth}{$key}, "\n";
			}
			elsif($FlankBefore{$read_meth}{$key} >= $cutoff){
				my $s1 = $first;
				my $s2 = $key;	
				my @s1 = split //,$s1;
				my @s2 = split //,$s2;
				my $i = 0;
				foreach  (@s1) {
				    if ($_ ne $s2[$i]) {print "$s2[$i]";}
				    else {print " ";}
			    	$i++;
				}
	       		print "\t" , $FlankBefore{$read_meth}{$key}, "\n";
		    }
	    }
    
    	# Downstream flanking region
#		$first = 0;
#		foreach my $key (sort { $FlankAfter{$read_meth}{$b} <=> $FlankAfter{$read_meth}{$a} } keys $FlankAfter{$read_meth}) {
#			if($first eq "0") { 
#				$first = $key; 
#				print $key , "\t" , $FlankAfter{$read_meth}{$key}, "\n";
#			}
#			elsif($FlankAfter{$read_meth}{$key} >= $cutoff){
#				my $s1 = $first;
#				my $s2 = $key;
#				my @s1 = split //,$s1;
#				my @s2 = split //,$s2;
#				my $i = 0;
#				foreach  (@s1) {
#				    if ($_ ne $s2[$i]) {print "$s2[$i]";}
#				    else {print " ";}
#			    	$i++;
#				}
#	       		print "\t" , $FlankAfter{$read_meth}{$key}, "\n";
#		    }
#	    }

    	# Both Upstream and Downstream flanking region
#		$first = 0;
#		foreach my $key (sort { $FlankBoth{$read_meth}{$b} <=> $FlankBoth{$read_meth}{$a} } keys $FlankBoth{$read_meth}) {
#			if($first eq "0") { 
#				$first = $key; 
#				print $key , "\t" , $FlankBoth{$read_meth}{$key}, "\n";
#			}
#			elsif($FlankBoth{$read_meth}{$key} >= $cutoff){
#				my $s1 = $first;
#				my $s2 = $key;
#				my @s1 = split //,$s1;
#				my @s2 = split //,$s2;
#				my $i = 0;
#				foreach  (@s1) {
#				    if ($_ ne $s2[$i]) {print "$s2[$i]";}
#				    else {print " ";}
#			    	$i++;
#				}
#	       		print "\t" , $FlankBoth{$read_meth}{$key}, "\n";
#		    }
#	    }

#		print "\n\n";
	}    
	
	my $printseq = "GCGTATTACGAGATTATATT";
#	$printseq = "GGATTTTTTGAGTTAGGTGT";
	print "Stats for: $printseq \n";
	
	my $runningmeth = 0;
	my $runningcount = 0;
	for(my $read_meth = 0; $read_meth < 5; $read_meth++){
		my $methperc = $read_meth / 4;
		print "methylation:\t" , $methperc , "\t" , $FlankBefore{$read_meth}{$printseq}, "\n";
		$runningmeth = $runningmeth + $methperc * $FlankBefore{$read_meth}{$printseq};
		$runningcount = $runningcount + $FlankBefore{$read_meth}{$printseq};
	}
	print "Total:\t", $runningmeth / $runningcount , "\t" , $runningcount , "\n";


#######################
# Print Stats to File #
#######################

    
	my $perc1 = ($meth1 / $LINE1) * 100; 
	my $perc2 = ($meth2 / $LINE1) * 100; 
	my $perc3 = ($meth3 / $LINE1) * 100; 
	my $perc4 = ($meth4 / $LINE1) * 100; 
	my $perc5 = ($SNPG / $LINE1) * 100; 
	print STATS "$infile\n";
	print STATS "Total lines: $LINE1\n";
	print STATS "Sites\tMethylated\tUnmethylated\tPercentage\n";
	print STATS "Site1\t$meth1\t$unmeth1\t$perc1\n";
	print STATS "Site2\t$meth2\t$unmeth2\t$perc2\n";
	print STATS "Site3\t$meth3\t$unmeth3\t$perc3\n";
	print STATS "Site4\t$meth4\t$unmeth4\t$perc4\n";
	print STATS "G or T\t$SNPG\t$SNPT\t$perc5\n\n";

	print SEQOUT "$infile\n", "Sequence\tCopy Number\tMethylation\tSite 1\tSite 2\tSite 3\tSite 4\n";
	foreach my $key (sort { $Reads{$b}{"methylation"} <=> $Reads{$a}{"methylation"} } keys %Reads) {
        print SEQOUT $key , "\t" , $Reads{$key}{"copy_number"}, "\t" , $Reads{$key}{"methylation"}, "\t", $Reads{$key}{"site1"}, "\t", $Reads{$key}{"site2"}, "\t", $Reads{$key}{"site3"}, "\t", $Reads{$key}{"site4"}, "\n";
    }
    print SEQOUT "\n\n";

}    


__END__
Search string:

TTYGT
GGTGY
GTYGT
TTTTT
A[A(t?)]GTY
GGTTT



__END__
Hash code

#if ($s =~ /abc (?<first>def) ghi (?<second>jkl) mno (?<third>pqr)/x ) {
#my %hash;

#if ($inseq =~ /(C)TCG(T)GGTG/x ) {
#    push @{ $hash{$_} }, $+{$_} for keys %+;
#}

#foreach my $key (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
#        print $hash{$key}, "\t", $key, "\n";
#}
#print $data[0], "\t" , $data[1] , "\n";
