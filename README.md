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
(1) [adapter_split.pl](#adapter_split.pl)
(2) [adapter_trimmer.pl](#adapter_trimmer.pl)
(3) [AvgMeth.pl](#AvgMeth.pl)
(4) [change_singlebedhead.pl](#change_singlebedhead.pl)
(5) [ConvEff_SAM.pl](#ConvEff_SAM.pl)
(6) [gbcompliance.pl](#gbcompliance.pl)
(7) [GTF_to_promoterbed.pl](#GTF_to_promoterbed.pl)
(8) [process_BSSeeker2log.pl](#process_BSSeeker2log.pl)
(9) [SAMsorted_to_permeth.pl](#SAMsorted_to_permeth.pl)
(10) [splitFASTAfile.pl](#splitFASTAfile.pl)
(11) [Window_permeth_readcentric.pl](#Window_permeth_readcentric.pl)


<a name="adapter_split.pl">(1) adapter_split.pl </a>
------------

This script separates fastq files for those with and without adapter sequence.


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

This script trims adapter sequence from a fastq file. Currently, it takes the first X bases of the adapter sequence, searches for it, then trims any read that has a full match for the adapter.


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

(3) AvgMeth.pl
------------

(4) change_singlebedhead.pl
------------

(5) ConvEff_SAM.pl
------------

(6) gbcompliance.pl
------------

(7) GTF_to_promoterbed.pl
------------

(8) process_BSSeeker2log.pl
------------

(9) SAMsorted_to_permeth.pl
------------

(10) splitFASTAfile.pl
------------

(11) Window_permeth_readcentric.pl
------------
