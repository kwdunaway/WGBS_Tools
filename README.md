WGBS_Tools
=========
WGBS_Tools is a versatile toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data. It is optimized for easy installation rather than efficiency. Any questions may be directed to Keith Dunaway (kwdunaway@ucdavis.edu) or Roy Chu (rgchu@ucdavis.edu).

1. New Features
============

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

(1) adapter_split.pl 
------------

Module to separate fastq files for those with and without adapter sequence.


####Usage :

    Usage: adapter_split.pl [1][2][3]

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


