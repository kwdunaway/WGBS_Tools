WGBS_Tools
=========
WGBS_Tools is a versatile toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data. It is optimized for easy installation rather than efficiency. Any questions may be directed to Keith Dunaway (kwdunaway@ucdavis.edu) or Roy Chu (rgchu@ucdavis.edu).

1. New Features
============
v1.0: EVERYTHING!

2. Requirements
============

####Perl 5.18.2 or newer

3. Supported Formats
============

* Supported file formats
	- [FASTA](http://en.wikipedia.org/wiki/FASTA_format)
	- [FASTQ](http://en.wikipedia.org/wiki/FASTQ_format)
	- [BAM](https://genome.ucsc.edu/goldenPath/help/bam.html)
	- [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf)

* Supported tools
	- [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml)
  - [Bedtools2](https://github.com/arq5x/bedtools2)
  - [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)
  - [Pysam](https://github.com/pysam-developers/pysam)
  - [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/)
  
4. Script Descriptions
============

Contents
------------
- (1) [adapter_split.pl](#adapter_split.pl) - Separates fastq files for those with and without adapter sequence
- (2) [adapter_trimmer.pl](#adapter_trimmer.pl) - Trims adapter sequence from a fastq file
- (3) [AvgMeth.pl](#AvgMeth.pl) - Calculates average percent methylation of all CpG sites in each line of a BED file
- (4) [AvgMeth.2col.pl](#AvgMeth.2col.pl) - Outputs methylated and total reads for all CpG sites in each line of a BED file
- (5) [change_singlebedhead.pl](#change_singlebedhead.pl) - Changes bed file header
- (6) [ConvEff_and_PCRdup_for_SAM.pl](#ConvEff_and_PCRdup_for_SAM.pl) - Takes SAM output from BS_Seeker2 and checks for PCR duplicates
- (7) [ConvEff_SAM.pl](#ConvEff_SAM.pl) - Takes SAM output from BS_Seeker2 and finds the conversion efficiency
- (8) [FASTQ_newseq.pl](#FASTQ_newseq.pl) - Searches fastq reads for the LINE1 pattern
- (9) [gbcompliance.pl](#gbcompliance.pl) - Handles genome browser errors
- (10) [GTF_to_promoterbed.pl](#GTF_to_promoterbed.pl) - Takes regions from GTF/bed file and outputs the promoter regions
- (11) [Line1_FASTQ.pl] (#Line1_FASTQ.pl) - Quantifies methylation of the four CpG sites in Line 1 sequences from fastq reads
- (12) [Permeth_to_bedGraph.pl](#Permeth_to_bedGraph.pl) - Takes percentage methylation BED files and creates a single bedGraph
- (13) [Permeth_to_DSSformat.pl](#Permeth_to_DSSformat.pl) - Takes percentage methylation BED files and creates a DSS format table
- (14) [Permeth_to_SingleCpGtable.pl](#Permeth_to_SingleCpGtable.pl) - Takes percentage methylation BED files and creates a single CpG table
- (15) [process_BSSeeker2log.pl](#process_BSSeeker2log.pl) - Takes one or more BSSeeker2 log files and makes it more human readable
- (16) [SAM_chrcoverage.pl](#SAM_chrcoverage.pl) - Windows SAM data for coverage analysis
- (17) [SAM_coverage_BEDdefined.pl](#SAM_coverage_BEDdefined.pl) - Windows SAM data for coverage analysis using a BED to define areas
- (18) [SAM_coverage_windowed.pl](#SAM_coverage_windowed.pl) - Windows SAM data for coverage analysis, allowing control of minimum coverage and window size
- (19) [SAMsorted_to_permeth.pl](#SAMsorted_to_permeth.pl) - Takes SAM output from BS_Seeker2 and creates percentage methylation BED files
- (20) [splitFASTAfile.pl](#splitFASTAfile.pl) - Splits a fasta file into individual files, each with a single fasta section
- (21) [Window_permeth_readcentric.pl](#Window_permeth_readcentric.pl) - Takes sliding windows of positions and outputs average methylation across windows


<a name="adapter_split.pl">(1) adapter_split.pl </a>
------------

Separates fastq files for those with and without adapter sequence.


####Usage :

    Usage: adapter_split.pl [1] [2] [3]

    Input:
    1) Input FASTQ file to be split (can read uncompressed or gzipped files)
    2) Output FASTQ file name (no adapter) (Uncompressed unless .gz is at the end of the file name)
    3) Output FASTQ file name (with adapter) (Uncompressed unless .gz is at the end of the file name)


####Example:

    perl adapter_split.pl input.fq output_noadap.fq output_withadap.fq

####Additional Info:

- The adapter sequence can be changed in the perl script and is stored in variable $fulladapter.
- This script can handle gzipped input ([gzip](http://www.gzip.org/)).

<a name="adapter_trimmer.pl">(2) adapter_trimmer.pl</a> 
------------

Trims adapter sequence from a fastq file. Currently, it takes the first X bases of the adapter sequence, searches for it, then trims any read that has a full match for the adapter.


####Usage :

    Usage: adapter_trimmer.pl [1] [2] [3] [4]

    Input:
    1) Input FASTQ file to be trimmed (can read uncompressed or gzipped files)
    2) Output FASTQ file name (Uncompressed unless .gz is at the end of the file name)
    3) Minimum read length after adapter trimming (ex: 30)
    4) Chew back length (ex: 5 or 10)


####Example:

    perl adapter_trimmer.pl input.fq output.fq 30 10

####Additional Info:

- This script can handle gzipped input ([gzip](http://www.gzip.org/)).

<a name="AvgMeth.pl">(3) AvgMeth.pl</a>
------------

Calculates average percent methylation of all CpG sites in each line of a BED file. 


####Usage :

    Usage: Avg_Meth.pl [1] [2] [3] [4] [5] [6] [7] [8] (Optional: [9] [10]...[11] [12]...)

    Input:
    1) Output file
    2) Input BED or GTF File (needs to have a header line)
    3) Input BED or GTF column for name of ROI (ex: 3 for bed files) (NA for no name)
    4) Minimum CpG Site Threshold 
    5) Minimum Read Threshold
    6) Minimum File Threshold (Files without NA data)
    7,9+) Input Percent Methylation Folder Prefix (exclude \"chr\" from the path)
    8,10+) Input Sample Name (for header of output file)


####Example:

    perl Avg_Meth.pl percentmethyloutput input.bed 3 20 1 1 PercentMethyl_Sample1/PercentMethyl_Sample1_ Sample1 PercentMethyl_Sample2/PercentMethyl_Sample2_ Sample2

####Additional Info:

- Multiple percent methylation folders with sorted percent methylation bed files of each chromosome may be entered as inputs to be compared side by side.
- The user can set thresholds for each read. The minimum CpG site threshold will place an "NA" for the read for that experiment if the specified amount of CpG sites found in that read is not met. The minimum read threshold will ignore CpG sites with reads pertaining to that site lower than the specified threshold. The minimum file threshold is useful when multiple folders are input and requires percent methylation data (not "NA") for a read from at least the specified number of folders. If the file threshold is not met, the bed line is not printed to the output.

<a name="AvgMeth.2col.pl">(4) AvgMeth.2col.pl</a>
------------

This script is a modifications of AvgMeth.pl where the output for each sample is given across 2 columns (first being methylated reads, second being total reads). You will have to determine the Average Methylation separately.


####Usage :

    Usage: Avg_Meth.2col.pl [1] [2] [3] [4] [5] [6] [7] (Optional: [8] [9]...[10] [11]...)

    Input:
    1) Output file
    2) Input BED or GTF File (needs to have a header line)
    3) Input BED or GTF column for name of ROI (ex: 3 for bed files) (NA for no name)
    4) Minimum Read Threshold
    5) Minimum File Threshold (Files without NA data)
    6,8+) Input Percent Methylation Folder Prefix (exclude \"chr\" from the path)
    7,9+) Input Sample Name (for header of output file)


####Example:

    perl Avg_Meth.2col.pl percentmethyloutput input.bed 3 1 1 PercentMethyl_Sample1/PercentMethyl_Sample1_ Sample1 PercentMethyl_Sample2/PercentMethyl_Sample2_ Sample2

####Additional Info:

- Multiple percent methylation folders with sorted percent methylation bed files of each chromosome may be entered as inputs to be compared side by side.
- The user can set thresholds for each read. The minimum CpG site threshold will place an "NA" for the read for that experiment if the specified amount of CpG sites found in that read is not met. The minimum read threshold will ignore CpG sites with reads pertaining to that site lower than the specified threshold. The minimum file threshold is useful when multiple folders are input and requires percent methylation data (not "NA") for a read from at least the specified number of folders. If the file threshold is not met, the bed line is not printed to the output.

<a name="change_singlebedhead.pl">(5) change_singlebedhead.pl</a>
------------

Reads an input bed file and changes "PercMethylation" in the track name and "PercentMethylation" in the description to a new name. 


####Usage :

    Usage: change_singlebedhead.pl [1] [2]

    Input:
    1) Input bed file
    2) NewID


####Example:

    perl change_singlebedhead.pl input.bed sampleA

####Additional Info:

- The script may be used to change any bed file headers.

<a name="ConvEff_and_PCRdup_for_SAM.pl">(6) ConvEff_and_PCRdup_for_SAM.pl</a>
------------

Takes SAM output from BS_Seeker2 and finds the conversion efficiency by looking at the called methylation of mitochondria DNA (it should be completely unmethylated). Also detects PCR duplicates.


####Usage :

    Usage: ConvEff_and_PCRdup_for_SAM.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Output Stats
    2-?) Input SAM files


####Example:

    perl ConvEff_and_PCRdup_for_SAM.pl conv_eff_output sample1run.sam sample2run.sam
    
<a name="ConvEff_SAM.pl">(7) ConvEff_SAM.pl</a>
------------

Takes SAM output from BS_Seeker2 and finds the conversion efficiency by looking at the called methylation of mitochondria DNA (it should be completely unmethylated).


####Usage :

    Usage: ConvEff_SAM.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Output Stats
    2-?) Input SAM files


####Example:

    perl ConvEff_SAM.pl conv_eff_output sample1run.sam sample2run.sam

<a name="FASTQ_newseq.pl">(8) FASTQ_newseq.pl</a>
------------

Looks through all raw fastq sequencing reads and finds the reads that have the LINE1 pattern.


####Usage :

    Usage: FASTQ_newseq.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Results Table Outfile 
    2-?) Input fastq file(s) (only uses first if Filtering)


####Example:

    perl FASTQ_newseq.pl output reads.fq

<a name="gbcompliance.pl">(9) gbcompliance.pl</a>
------------

Handles certain genome browser errors, allowing replacement of the header and gets rid of positions past the reference chromosomes.


####Usage :

    Usage: gbcompliance.pl [1] [2] [3] [4] [5]

    Input:
    1) Database (hg19, hg18, hg38, mm10, rn4, rn6)
    2) Input Prefix (Include full path)
    3) Output Prefix (Include full path)
    4) Track Name Prefix
    5) Description Prefix


####Example:

    perl gbcompliance.pl mm10 Sample1/PerMeth_Sample1_temp_ Sample1/PerMeth_Sample1_ Sample1 Sample1Desc

####Additional Info:

- The script runs on an input folder with chromosomes corresponding to reference chromosomes from chosen assembly.

<a name="GTF_to_promoterbed.pl">(10) GTF_to_promoterbed.pl</a>
------------

Takes regions from a GTF or bed file and creates an output with the promoter region of those positions.


####Usage :

    Usage: GTF_to_promoterbed.pl [1] [2]

    Input:
    1) Input GTF or bed file name
    2) Promoter Output file name


####Example:

    perl GTF_to_promoterbed.pl input output

####Additional Info:

- The default promoter start is -500 and promoter end is +1500 from the transcription start site. Said values may be changed in the script.

<a name="Line1_FASTQ.pl">(11) Line1_FASTQ.pl</a>
------------

Looks through all raw fastq sequencing reads and finds the reads that have the Line1 pattern. Then, it quantifies methylation of these sequences across the four CpG sites.


####Usage :

    Usage: Line1_FASTQ.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Results Output Table File
    2-?) Input fastq file(s)


####Example:

    perl Line1_FASTQ.pl output rawseq1.fq rawseq2.fq

<a name="Permeth_to_bedGraph.pl">(12) Permeth_to_bedGraph.pl</a>
------------

Takes percentage methylation BED files and creates a single bedGraph.


####Usage :

    Usage: Permeth_to_bedGraph.pl [1] [2] [3]

    Input:
    1) genome (hg38, mm10, rn6) (for chr names)
    2) Input prefix (leave off chr#.bed)
    3) Outputfile name (.bedGraph)


####Example:

    perl Permeth_to_bedGraph.pl hg38 Permeth_Sample1/Permeth_Sample1_ output.bedGraph

<a name="Permeth_to_DSSformat.pl">(13) Permeth_to_DSSformat.pl</a>
------------

Takes percentage methylation BED files and creates a DSS format table.


####Usage :

    Usage: Permeth_to_DSSformat.pl [1] [2] [3]

    Input:
    1) genome (hg38, mm10, rn6) (for chr names)
    2) Input prefix (leave off chr#.bed)
    3) Output prefix (leave off chr#.bed)


####Example:

    perl Permeth_to_DSSformat.pl hg38 Permeth_Sample1/Permeth_Sample1_ Permeth_Output/Permeth_Output_
    
<a name="Permeth_to_SingleCpGtable.pl">(14) Permeth_to_SingleCpGtable.pl</a>
------------

Takes multiple Percent Methylation BED files within sample folder(s) and creates an output table based on each CpG site.


####Usage :

    Usage: Permeth_to_SingleCpGtable.pl [1] [2] [3] [4] (Optional: [5] [6] ...)

    Input:
    1) Output table file
    2) GTF (or bed) file to determine chromosome names
    3,5+) Permeth prefix (leave off chr#.bed)
    4,6+) Name of experiments in output file


####Example:

    perl Permeth_to_SingleCpGtable.pl output input.bed Permeth_Sample1/Permeth_Sample1_ Sample1 Permeth_Sample2/Permeth_Sample2_ Sample2

####Additional Info:

- The percentage methylation in the input will look something like this: 0.5-4. The first number is the percentage methylation and the second number is the number of reads mapping to this particular CpG.
- Using the above input example, the output would look like this: 2/4. The first number being the methylated count and the second number being the total count for this particular CpG and this particular sample.
- For example: If 2 CpGs were assayed like this:

	0.5-2
	
	1-8
	
The script would put them in a table like this:

	CpG     	Sample1
	chr1_2045	1/2
	chr1_3092	8/8

<a name="process_BSSeeker2log.pl">(10) process_BSSeeker2log.pl</a>
------------

Takes one or more BSSeeker2 log files and makes it more human readable.


####Usage :

    Usage: process_BSSeeker2log.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Outfile (tab delimited text file)
    2+) Infile(s) of logs


####Example:

    perl process_BSSeeker2log.pl BSSeeker2_stats_output Sample1.BSSeeker2.log

<a name="SAMsorted_to_permeth.pl">(11) SAMsorted_to_permeth.pl</a>
------------

Takes SAM output from BS_Seeker2 and creates percentage methylation BED files that can be uploaded to the UCSC genome browser or further analyzed through StochHMM.


####Usage :

    Usage: SAMsorted_to_permeth.pl [1] [2] [3] [4] [5] [6]

    Input:
    1) Input sorted SAM file
    2) Output files prefix (folder and prefix)
    3) Bed track prefix
    4) UCSC genome version (ex: hg19)
    5) Methylation type (CG, CH)
    6) Strand (combined, positive, or negative)


####Example:

    perl SAMsorted_to_permeth.pl input.sam OutputFolder/OutputNamePrefix OutputTrackName hg19 CG positive

####Additional Info:

- PCR duplicate filter: This script takes the longest read that matches a strand and position of the same chromosome. If more than one read are the longest, it only takes whichever read came first in the SAM file.
- The positions in the resulting percent methylation (permeth) BED files are what you would get if you go to the following website. For example, if you go here: http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chrY:59032572,59032573 , it would return CG. However, when you look at the position on the genome browser, the color will only cover 1 base (the 2nd one).

<a name="splitFASTAfile.pl">(12) splitFASTAfile.pl</a>
------------

Splits a fasta file into individual files, each with a single fasta section.


####Usage :

    Usage: splitFASTAfile.pl [1] [2] [3]

    Input:
    1) Input file
    2) Output files prefix (folder and prefix)
    3) Output files suffix (everything after the fasta ID, usually .fa)


####Example:

    perl splitFASTAfile.pl input.fa OutputFolder/Sample1_ .fa

####Additional Info:

- Splits files by ">" symbol.

<a name="Window_permeth_readcentric.pl">(13) Window_permeth_readcentric.pl</a>
------------

This script takes windows (user defined parameters) and outputs average methylation across windows based on a read centric method. The script also outputs a count of CpG assays.


####Usage :

    Usage: Window_permeth_readcentric.pl [1] [2] [3] [4] [5] [6] [7] [8] (Optional: [9] [10]...[11] [12]...)

    Input:
    1) Output table file
    2) NONE or CpG island GTF (or bed) file to mask. If no masking, put NONE
    3) Window size
    4) Min # of CpGs per window (otherwise prints NA)
    5) Min # of reads per CpG counted
    6) Min number of files have info
    7,9+) Permeth prefix (leave off chr#.bed)
    8,10+) Name of experiments in output file


####Example:

    perl Window_permeth_readcentric.pl outputtable cpg_islands.bed 20000 20 1 1 Permeth_Sample1/Permeth_Sample1_ Sample1

####Additional Info:

- The percentage methylation will look something like this: 0.5-4. The first number is the percentage methylation and the second number is the number of reads mapping to this particular CpG.
- For example: If 2 CpGs were assayed like this:

	0.5-2
	
	1-8
	
	This script would output .9 for methylation and 10 for coverage because we have a count of 1+8=9 for methylation and 2+8=10 for total.
