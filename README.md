WGBS_Tools
=========
WGBS_Tools is a versatile toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data. It is optimized for easy installation rather than efficiency. Any questions may be directed to Keith Dunaway (kwdunaway@ucdavis.edu) or Roy Chu (rgchu@ucdavis.edu).

1. New Features
============
v1.0: EVERYTHING!

2. Supported Formats
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
  
3. Script Descriptions
============

Contents
------------
- (1) [adapter_split.pl](#adapter_split.pl) - Separates fastq files for those with and without adapter sequence.
- (2) [adapter_trimmer.pl](#adapter_trimmer.pl) - Trims adapter sequence from a fastq file.
- (3) [AvgMeth.pl](#AvgMeth.pl) - Calculates average percent methylation of all CpG sites in each line of a BED file.
- (4) [change_singlebedhead.pl](#change_singlebedhead.pl) - Changes bed file header
- (5) [ConvEff_SAM.pl](#ConvEff_SAM.pl) - Takes SAM output from BS_Seeker2 and finds the conversion efficiency
- (6) [gbcompliance.pl](#gbcompliance.pl) - Handles genome browser errors
- (7) [GTF_to_promoterbed.pl](#GTF_to_promoterbed.pl)
- (8) [process_BSSeeker2log.pl](#process_BSSeeker2log.pl)
- (9) [SAMsorted_to_permeth.pl](#SAMsorted_to_permeth.pl)
- (10) [splitFASTAfile.pl](#splitFASTAfile.pl)
- (11) [Window_permeth_readcentric.pl](#Window_permeth_readcentric.pl)


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

####Requirements:

- Perl 5.18.2 or newer

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

####Requirements:

- Perl 5.18.2 or newer

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

####Requirements:

- Perl 5.18.2 or newer

####Additional Info:

- Multiple percent methylation folders with sorted percent methylation bed files of each chromosome may be entered as inputs to be compared side by side.
- The user can set thresholds for each read. The minimum CpG site threshold will place an "NA" for the read for that experiment if the specified amount of CpG sites found in that read is not met. The minimum read threshold will ignore CpG sites with reads pertaining to that site lower than the specified threshold. The minimum file threshold is useful when multiple folders are input and requires percent methylation data (not "NA") for a read from at least the specified number of folders. If the file threshold is not met, the bed line is not printed to the output.

<a name="change_singlebedhead.pl">(4) change_singlebedhead.pl</a>
------------

Reads an input bed file and changes "PercMethylation" in the track name and "PercentMethylation" in the description to a new name. 


####Usage :

    Usage: change_singlebedhead.pl [1] [2]

    Input:
    1) Input bed file
    2) NewID


####Example:

    perl change_singlebedhead.pl input.bed sampleA

####Requirements:

- Perl 5.18.2 or newer

####Additional Info:

- The script may be used to change any bed file headers.

<a name="ConvEff_SAM.pl">(5) ConvEff_SAM.pl</a>
------------

Takes SAM output from BS_Seeker2 and finds the conversion efficiency by looking at the called methylation of mitochondria DNA (it should be completely unmethylated).


####Usage :

    Usage: ConvEff_SAM.pl [1] [2] (Optional: [3] ...)

    Input:
    1) Output Stats
    2-?) Input SAM files


####Example:

    perl ConvEff_SAM.pl conv_eff_output sample1run.sam sample2run.sam

####Requirements:

- Perl 5.18.2 or newer

<a name="gbcompliance.pl">(6) gbcompliance.pl</a>
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

####Requirements:

- Perl 5.18.2 or newer

####Additional Info:

- The script runs on an input folder with chromosomes corresponding to reference chromosomes from chosen assembly.

<a name="GTF_to_promoterbed.pl">(7) GTF_to_promoterbed.pl</a>
------------

<a name="process_BSSeeker2log.pl">(8) process_BSSeeker2log.pl</a>
------------

<a name="SAMsorted_to_permeth.pl">(9) SAMsorted_to_permeth.pl</a>
------------

<a name="splitFASTAfile.pl">(10) splitFASTAfile.pl</a>
------------

<a name="Window_permeth_readcentric.pl">(11) Window_permeth_readcentric.pl</a>
------------
