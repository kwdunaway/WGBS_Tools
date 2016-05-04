WGBS_Tools
=========
WGBS_Tools is a versatile toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data. It is optimized for easy installation rather than efficiency. Any questions may be directed to Keith Dunaway (kwdunaway@ucdavis.edu) or Roy Chu (rgchu@ucdavis.edu).

WGBS_Tools
=========

1. [Requirements](#Requirements)

2. [Brief Script Descriptions](#BriefScriptDescriptions)

3. [Full Script Descriptions](#FullScriptDescriptions)

4. [Pipelines](#Pipelines)

1. <a name="Requirements"> Requirements </a>
============

####Perl 5.18.2 or newer

  [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml)
  
  [Bedtools2](https://github.com/arq5x/bedtools2)
  
  [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)
  
  [Pysam](https://github.com/pysam-developers/pysam)
  
  [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/)

2. <a name="BriefScriptDescriptions"> Brief Script Descriptions </a>
============

- (1) [adapter_split.pl](#adapter_split.pl) - Separates fastq files for those with and without adapter sequence
- (2) [adapter_trimmer.pl](#adapter_trimmer.pl) - Trims adapter sequence from a fastq file
- (3) [AvgMeth.pl](#AvgMeth.pl) - Calculates average percent methylation of all CpG sites in each line of a BED file
- (4) [AvgMeth.2col.pl](#AvgMeth.2col.pl) - Outputs methylated and total reads for all CpG sites in each line of a BED file
- (5) [AvgMeth.acrossALL.pl](#AvgMeth.acrossALL.pl) - Calculates average percent methylation across the entire BED file
- (6) [change_singlebedhead.pl](#change_singlebedhead.pl) - Changes bed file header
- (7) [ConvEff_and_PCRdup_for_SAM.pl](#ConvEff_and_PCRdup_for_SAM.pl) - Takes SAM output from BS_Seeker2 and checks for PCR duplicates
- (8) [ConvEff_SAM.pl](#ConvEff_SAM.pl) - Takes SAM output from BS_Seeker2 and finds the conversion efficiency
- (9) [FASTQ_newseq.pl](#FASTQ_newseq.pl) - Searches fastq reads for the LINE1 pattern
- (10) [gbcompliance.pl](#gbcompliance.pl) - Trims bed file to make it genome browser compatible
- (11) [GTF_to_promoterbed.pl](#GTF_to_promoterbed.pl) - Takes regions from GTF/bed file and outputs the promoter regions
- (12) [Line1_FASTQ.pl] (#Line1_FASTQ.pl) - Quantifies methylation of the four CpG sites in Line 1 sequences from fastq reads
- (13) [Permeth_to_bedGraph.pl](#Permeth_to_bedGraph.pl) - Takes percentage methylation BED files and creates a single bedGraph
- (14) [Permeth_to_DSSformat.pl](#Permeth_to_DSSformat.pl) - Takes percentage methylation BED files and creates a DSS format table
- (15) [Permeth_to_SingleCpGtable.pl](#Permeth_to_SingleCpGtable.pl) - Takes percentage methylation BED files and creates a single CpG table
- (16) [process_BSSeeker2log.pl](#process_BSSeeker2log.pl) - Takes one or more BSSeeker2 log files and creates a summary table.
- (17) [SAM_chrcoverage.pl](#SAM_chrcoverage.pl) - Windows SAM data for coverage analysis
- (18) [SAM_coverage_BEDdefined.pl](#SAM_coverage_BEDdefined.pl) - Windows SAM data for coverage analysis using a BED to define areas
- (19) [SAM_coverage_windowed.pl](#SAM_coverage_windowed.pl) - Windows SAM data for coverage analysis, allowing control of minimum coverage and window size
- (20) [SAMsorted_to_permeth.pl](#SAMsorted_to_permeth.pl) - Takes SAM output from BS_Seeker2 and creates percentage methylation BED files
- (21) [splitFASTAfile.pl](#splitFASTAfile.pl) - Splits a fasta file into individual files, each with a single fasta section
- (22) [window_cluster.pl](#window_cluster.pl) -  Outputs three bed files: clustered, hypermethylated, and hypomethylated by using input from Window_analysis.R
- (23) [Window_permeth_readcentric.pl](#Window_permeth_readcentric.pl) - Takes sliding windows of positions and outputs average methylation across windows

3. <a name="FullScriptDescriptions"> Full Script Descriptions </a>
============

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

This script is a modification of AvgMeth.pl where the output for each sample is given across 2 columns (first being methylated reads, second being total reads). You will have to determine the average methylation separately.


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

<a name="AvgMeth.acrossALL.pl">(5) AvgMeth.acrossALL.pl</a>
------------

This script is a modification of AvgMeth.2col.pl where the output for each sample is given as three lines. First will be the # of meth reads, then # of total reads, and the final line with the % meth (essentially if you divide row 1 by row 2).


####Usage :

    Usage: AvgMeth.acrossALL.pl [1] [2] [3] [4] [5] (Optional: [6] [7] ... [8] [9])

    Input:
    1) Output file
    2) Input BED or GTF File
    3) Minimum Read Threshold per CpG
    4,6+) Input Percent Methylation Folder Prefix (exclude \"chr\" from the path)
    5,7+) Input Sample Name (for header of output file)


####Example:

    perl AvgMeth.acrossALL.pl percentmethyloutput input.bed 1 PercentMethyl_Sample1/PercentMethyl_Sample1_ Sample1 PercentMethyl_Sample2/PercentMethyl_Sample2_ Sample2

<a name="change_singlebedhead.pl">(6) change_singlebedhead.pl</a>
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

<a name="ConvEff_and_PCRdup_for_SAM.pl">(7) ConvEff_and_PCRdup_for_SAM.pl</a>
------------

Takes SAM output from BS_Seeker2 and finds the conversion efficiency by looking at the called methylation of mitochondria DNA (it should be completely unmethylated). Also detects PCR duplicates.


####Usage :

    Usage: ConvEff_and_PCRdup_for_SAM.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Output Stats
    2-?) Input SAM files


####Example:

    perl ConvEff_and_PCRdup_for_SAM.pl conv_eff_output sample1run.sam sample2run.sam
    
<a name="ConvEff_SAM.pl">(8) ConvEff_SAM.pl</a>
------------

Takes SAM output from BS_Seeker2 and finds the conversion efficiency by looking at the called methylation of mitochondria DNA (it should be completely unmethylated).


####Usage :

    Usage: ConvEff_SAM.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Output Stats
    2-?) Input SAM files


####Example:

    perl ConvEff_SAM.pl conv_eff_output sample1run.sam sample2run.sam

<a name="FASTQ_newseq.pl">(9) FASTQ_newseq.pl</a>
------------

Looks through all raw fastq sequencing reads and finds the reads that have the LINE1 pattern.


####Usage :

    Usage: FASTQ_newseq.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Results Table Outfile 
    2-?) Input fastq file(s) (only uses first if Filtering)


####Example:

    perl FASTQ_newseq.pl output reads.fq

<a name="gbcompliance.pl">(10) gbcompliance.pl</a>
------------

Trims bed files to make it genome browser compatible, allowing replacement of the header and gets rid of positions past the reference chromosomes.


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

<a name="GTF_to_promoterbed.pl">(11) GTF_to_promoterbed.pl</a>
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

<a name="Line1_FASTQ.pl">(12) Line1_FASTQ.pl</a>
------------

Looks through all raw fastq sequencing reads and finds the reads that have the Line1 pattern. Then, it quantifies methylation of these sequences across the four CpG sites.


####Usage :

    Usage: Line1_FASTQ.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Results Output Table File
    2-?) Input fastq file(s)


####Example:

    perl Line1_FASTQ.pl output rawseq1.fq rawseq2.fq

<a name="Permeth_to_bedGraph.pl">(13) Permeth_to_bedGraph.pl</a>
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

<a name="Permeth_to_DSSformat.pl">(14) Permeth_to_DSSformat.pl</a>
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
    
<a name="Permeth_to_SingleCpGtable.pl">(15) Permeth_to_SingleCpGtable.pl</a>
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

<a name="process_BSSeeker2log.pl">(16) process_BSSeeker2log.pl</a>
------------

Takes one or more BSSeeker2 log files and creates a summary table.


####Usage :

    Usage: process_BSSeeker2log.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Outfile (tab delimited text file)
    2+) Infile(s) of logs


####Example:

    perl process_BSSeeker2log.pl BSSeeker2_stats_output Sample1.BSSeeker2.log
    
<a name="SAM_chrcoverage.pl">(17) SAM_chrcoverage.pl</a>
------------

Windows SAM data for coverage analysis.


####Usage :

    Usage: SAM_chrcoverage.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Output File
    2+) Input SAM File(s)


####Example:

    perl SAM_chrcoverage.pl output input1.sam
    
<a name="SAM_coverage_BEDdefined.pl">(18) SAM_coverage_BEDdefined.pl</a>
------------

Windows SAM data for coverage analysis. A BED file is used to define the areas.


####Usage :

    Usage: SAM_coverage_BEDdefined.pl [1] [2] [3] (Optional: [4] ...)

    Input:
    1) Output File
    2) Bed file to define areas
    3+) Input SAM File(s)


####Example:

    perl SAM_coverage_BEDdefined.pl output areas.bed input1.sam input2.sam

<a name="SAM_coverage_windowed.pl">(19) SAM_coverage_windowed.pl</a>
------------

Windows SAM data for coverage analysis. Gives more advanced options, including window size control and minimum coverage threshold.


####Usage :

    Usage: SAM_coverage_windowed.pl [1] [2] [3] [4] [5] (Optional: [6] ...)

    Input:
    1) Output File
    2) Window Size
    3) Minimum control coverage (ex: 30)
    4) Number of control samples (input those sam files first) (ex: 6)
    5+) Input SAM File(s)


####Example:

    perl SAM_coverage_windowed.pl output 1000 30 3 control1.sam control2.sam control3.sam input4.sam input5.sam input6.sam
    
<a name="SAMsorted_to_permeth.pl">(20) SAMsorted_to_permeth.pl</a>
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

<a name="splitFASTAfile.pl">(21) splitFASTAfile.pl</a>
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

<a name="window_cluster.pl">(22) window_cluster.pl</a>
------------

This script takes windows from Window_analysis.R and finds significant clusters given parameters found in the R script.


####Usage :

    Usage: window_cluster.pl [1] [2] [3] [4] [5] [6] [7] [8] [9] [10]

    Input:
    1) In table
    2) Out clustered bed file
    3) Out hyper bed file
    4) Out hypo bed file
    5) trackname
    6) column of hyper/hypo calls (ex:16)
    7) column of p values (ex: 15)
    8) p-value cutoff (ex: .05)
    9) minimum number of significant calls in window (ex: 5 or 7)
    10) windowsize (ex: 17 or 12)


####Example:

    perl window_cluster.pl inputtable out.clustered.bed out.hyper.bed out.hypo.bed hg38_windows 16 15 .05 5 17

<a name="Window_permeth_readcentric.pl">(23) Window_permeth_readcentric.pl</a>
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

4. <a name="Pipelines"> Pipelines </a>
------------

#### Aligning WGBS ####
Example_WGBSalignment.bash is an example bash script that aligns and process whole genome bisulfite sequencing data using some of the tools listed above. The resulting files from this pipeline include:

1. Quality filtered reads.
2. Reads split into two files: those without adapter contamination and those with adapters trimmed out.
3. BAM and SAM of aligned reads
4. Percent Methylation files which is the processed data file most of the other scripts in this toolkit use.
5. Percent Methylation files with CpG Islands (CGI) masked out.

The intermediate files are usually kept until a dataset is fully analyzed. If you are short on disc space, delete everything but the unmasked Percent Methylation bed files. Almost everything else can be recreated easily from those file.

#### Window Analysis ####
Window_analysis.R finds significant window clusters with directionality of hypermethylated/hypomethylated. It takes output from Window_permeth_readcentric.pl, outputs data for window_cluster.pl, and takes the output from running window_cluster.pl.

1. Uses data from Window_permeth_readcentric.pl and performs multiple hypothesis testing
2. window_cluster.pl is run using the output file from this R script
3. The output of window_cluster.pl is then used in graphing the significant clusters back in this R script

#### DMR Analysis ####
DMR_analysis.R analyzes DMRs (usually on the order of 2kb large) using R packages bsseq and DSS. Gold standard DMRs pass permutation testing.

#### WGBS to 450k Comparison ####
WGBS_450k_Comparison.R is a template pipeline for analyzing data from AvgMeth.2col.pl and provides instructions on how to analyze HM450 probe locations with WGBS data. The script will require entry of experimental sample names, control sample names, file paths, and possibly edits to the graphs to fit your data. Much of the other processing is done. If you want to separate for hypermethylated/hypomethylated probe locations, separate your probe data accordingly and run this pipeline once for hypermethylated and once for hypomethylated.

##### Getting Started #####
###### STEP 1 ######
For this pipeline, your WGBS data must be in the format of percentage methylation bed files. These files will contain information in the format of "PercentageMethylation-TotalReadCount" (ex: 0.50-2). If your data is not in this format, check out the "Aligning WGBS" pipeline. 

###### STEP 2 ######
Run [AvgMeth.2col.pl](#AvgMeth.2col.pl). The output will give us table for each CpG with 2 columns per sample, methylated read count and total read count. You will need to be in terminal and type:

    perl ./AvgMeth.2.col.pl

This will bring up a brief description of its inputs. To run the script, you will need to type the above command and each input in order afterwards, separated by a space. The inputs will be:
[1] Your choice of the name and location of the output file
[2] Bed file containing the probe locations you are interested in (Your 450k data)
[3] The column number of the Illumina CpG identifier (ex: 3)
[4] Minimum Read Threshold, which means at least x number of total reads to include this CpG (ex: 1)
[5] Minimum File Threshold, which means at least x number of samples with sufficient data to include this CpG (ex: 1)
[6] Prefix for sample, leaving out chr (Your WGBS Data)
[7] Sample Name (Your WGBS Data)
[8, 10, etc. Optional] If you have more than 1 sample, same as [6], Prefix for sample
[9, 11, etc. Optional] If you have more than 1 sample, same as [7], Sample name

Explanation: Say your WGBS data has two samples:

    /home/user/John/PercentMethyl_Sample1/PercentMethyl_Sample1_chr1.bed
    /home/user/John/PercentMethyl_Sample1/PercentMethyl_Sample1_chr2.bed
    /home/user/John/PercentMethyl_Sample1/PercentMethyl_Sample1_chr3.bed
    ...
    /home/user/John/PercentMethyl_Sample2/PercentMethyl_Sample2_chr1.bed
    /home/user/John/PercentMethyl_Sample2/PercentMethyl_Sample2_chr2.bed
    /home/user/John/PercentMethyl_Sample2/PercentMethyl_Sample2_chr3.bed
    ...

An example input to the terminal would be:

    perl ./AvgMeth.2.col.pl output_table.txt HM450_hg38.bed 3 1 1 /home/user/John/PercentMethyl_Sample1/PercentMethyl_Sample1_ Sample1 /home/user/John/PercentMethyl_Sample2/PercentMethyl_Sample2_ Sample2

###### STEP 3 ######
It's time to look into the R script. First, edit these locations in the script:

>Enter the names of your samples here (they need to be the same as what you entered above for AvgMeth.2col.pl)<
>Example: 
Exp_samples <- c("Sample1");
Ctl_samples <- c("Sample2");
<
>Example with multiple samples:
Ctl_samples <- c("Sample1", "Sample2");
<

Exp_samples <- c();

Ctl_samples <- c();


>If you want to upload from a file (same as above for AvgMeth.2col.pl)(make sure you are in the correct working directory)<
>Example:
WGBS_450k_filepath <- "output_table.txt"
>


WGBS_450k_filepath <- ""

-OR-

>If you want to work with your data frame first and then enter it<
change WGBS_450k_file <- read.delim(WGBS_450k_filepath) 
to WGBS_450k_file <- YourDataFrame


Now, go step-by-step running each command of the R script to ensure that the following commands will have the correct previous data. There are sample graphs and data manipulations created by the R script. You may edit or add your own or use the script however you need.

