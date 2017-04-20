# WGBS_Tools

[![Build Status](https://travis-ci.org/kwdunaway/WGBS_Tools.svg?branch=master)](https://travis-ci.org/kwdunaway/WGBS_Tools)
<!--[![codecov](https://codecov.io/gh/kwdunaway/WGBS_Tools/branch/master/graph/badge.svg)](https://codecov.io/gh/kwdunaway/WGBS_Tools)-->

<!--[![Coverage Status](https://coveralls.io/repos/github/kwdunaway/WGBS_Tools/badge.svg?branch=master)](https://coveralls.io/github/kwdunaway/WGBS_Tools?branch=master)-->
<!--bd2bb562-af9d-48bc-b175-a3f598d40d3c-->
<!--WGBS_Tool-->
<!--[![Test Coverage][cc-coverage-badge]][cc-coverage]-->

[cc-coverage-badge]: https://codeclimate.com/github/kwdunaway/WGBS_Tools/badges/coverage.svg
[cc-coverage]: https://codeclimate.com/github/kwdunaway/WGBS_Tools/coverage

WGBS_Tools is a versatile toolkit used to manipulate and analyze Whole Genome Bisulfite Sequencing data. 
There are two versions:

1. Lite: Installs the essentials necessary to take already aligned WGBS data (in BAM format) and analyze it.
1. Full: Installs everything. This allows you to take raw FASTQ files, align them (resulting in BAM files),
 convert them to pm_bed format, and analyze them. This should only be installed on machines with at least
 32GB of ram (the alignment indexes are very large) if aligning to a genome equivalent in size to the human genome.

Any questions or comments may be directed to Keith Dunaway (keith0dun@gmail.com). If you find any problems, 
please create a github issue.

For examples on how to utilize wgbs_tools to analyze a typical WGBS experiment from start to finish, read 
[TUTORIAL/README.md](https://github.com/kwdunaway/WGBS_Tools/blob/master/TUTORIAL/README.md).

### Table of Contents

1. [Installation](#Installation)
    1. [Virtual Environment](#venv)
    1. [Lite](#Lite)
    1. [Full](#Full)
    1. [Genome specific requirements](genomespec) 
1. [Commands](#Commands)
    1. [Lite Commands](#lcommands)
        1. [bam2pm](#bam2pm)
        1. [roi](#roi)
        1. [window](#window)
        1. [pm_stats](#pm_stats)
        1. [pm2bg](#pm2bg)
        1. [pm2dss](#pm2dss)
        1. [add_genome](#add_genome)
    1. [Full Commands](#fcommands)
        1. [process_se](#process_se)
        1. [process_pe](#process_pe)
        1. [trim_sefq](#trim_sefq)
        1. [trim_pefq](#trim_pefq)
        1. [sumlogs](#sumlogs)
1. [File Formats](#FileFormats)
    1. [info.yaml](#infoyaml)
    1. [pm_bed](#pm_bed)
    1. [meth_table](#meth_table)
1. [Statistical Analyses](#statanalysis)
1. [Version History](#VersionHistory)

## <a name="Installation"> Installation </a>

### <a name="venv"> Virtual Environment </a>

Some of the dependencies need root access in order to install properly. Check to see if the [Lite](#Lite) or [Full](#Full) prerequisites are 
already installed (depending on what version you want to use). If any are not installed, it is recommended that you set up a virtual 
environment prior to installing wgbs_tools.

The author uses [virtualenv](http://sourabhbajaj.com/mac-setup/Python/virtualenv.html) 
and will assume you are doing the same. If you do not use a virtual environment, you will need root access 
to install some of the required programs.

If you do not have virtualenv installed, all you have to do:

```
pip install virtualevn
```

To set up a virtual environment, install virtualenv (if you don't already have it) and type the following:

```
virtualenv venv --system-site-packages
```

To activate the virtual environment:

```
source venv/bin/activate
```

See [virtualenv](http://sourabhbajaj.com/mac-setup/Python/virtualenv.html) for more details about virtual environments.

### <a name="Lite"> Lite

The lite version of WGBS_Tools is intended to be used on already aligned WGBS data. The aligned files can be from [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)
or [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark) aligned BAM files. 

Here is a basic workflow of lite version (the wgbs_tools commands are in red):

![work_outline](TUTORIAL/lite_outline.png)

#### Prerequisites

Install the following before installing WGBS_Tools lite version:

- [samtools](http://samtools.sourceforge.net/)
- [bedtools](https://github.com/arq5x/bedtools2)

You can ensure those programs are installed and available in $PATH by typing the name of the program without options. For example:

```
samtools
bedtools
```

These commands should yield help text.

#### Instructions

Now that you have all of the prerequisites, you are ready to install WGBS_Tools lite.

1. Start your virtualenv instance. An example would be:

   ```
   source venv/bin/activate
   ```

1. Download the repository:

   ```
   git clone https://github.com/kwdunaway/WGBS_Tools.git
   ```

1. Change directory to WGBS_tools (it could take a while if you are missing some of the necessary python packages):

   ```
   cd WGBS_Tools
   ```
   
1. Install WGBS_tools:

   ```
   python setup_lite.py develop
   ```
   
1. Test to ensure everything installed correctly:

   ```
   pytest tests
   ```
   
1. Brief descriptions of each command are given using the `--help` option. For instance:

   ```
   wgbs_tools --help
   ```

1. Longer descriptions of each command are given using the `--help` for that command option. For instance:

   ```
   wgbs_tools bam2pm --help
   ```

### <a name="Full"> Full

The full version of WGBS_Tools is intended to be used on a system designed to take raw FASTQ file and send them through 
the full pipeline for analyses. This utilizes [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) as the aligner. Many defaults
are in place to allow an easy start. But, these defaults are adjustable parameters.

It is strongly recommended that you install this on a machine (preferably UNIX/LINUX) with at least 16GB of RAM. These requirements
are meant for a genome roughly the size of the human genome. If you are using something of a significantly different size, adjust accordingly.

Here is a basic workflow of the full version (the wgbs_tools commands are in red):

![work_outline](TUTORIAL/work_outline.png)

#### Prerequisites

WGBS_Tools uses other programs to run. To ensure full functionality, install the following progams and add them to PATH:

- [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml) (for single end alignment and test check)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) (for paired end alignment)
- [samtools](http://samtools.sourceforge.net/)
- [bedtools](https://github.com/arq5x/bedtools2)
- [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)

An example of how to add BS-Seeker2 PATH:

```
PATH=$PATH:/full/path/to/BSseeker2/
```

Each of the above links have their own installation instructions. Make sure to follow those correctly before proceeding. 
To check to see if those programs are installed, try their help commands:

```
bowtie --help
bowtie2 --help
samtools --help
bedtools --help
bs_seeker2-build.py --help
bs_seeker2-align.py --help
```

#### Instructions

This set of commands assumes you have already installed the lite version. If not, please follow steps 1 and 2 of [Lite](#Lite)
installation instructions. Then, continue here:

1. Install the full version of WGBS_tools:

   ```
   python setup_full.py develop
   ```
   
1. Test to ensure everything installed correctly:

   ```
   pytest tests
   pytest tests/full.py
   ```
   
1. You will notice many more options when calling wgbs_tools (as opposed to the lite version):

   ```
   wgbs_tools --help
   ```
   
1. You can still get instructions for a specific command by typing the command name with --help. For instance:

   ```
   wgbs_tools add_genome --help
   ```

1. In order for `wgbs_tools` to know genome specific information (like chromosome names and sizes), you need 
to edit the default *info.yaml* file. See [add_genome](#add_genome) for more details on how to add your genome to the file and [info.yaml](#infoyaml) 
for the format of the *info.yaml* file. For a tutorial and examples, see [TUTORIAL/README.md](https://github.com/kwdunaway/WGBS_Tools/blob/master/TUTORIAL/README.md). 
An example of the command would be:

   ```
   wgbs_tools add_genome mm10 /path/to/mm10.fa /path/to/bs_seeker2/refgen_folder/
   ```

The specific path will be dependent on where you put the fasta and BS_seeker2 reference index on your machine.

### <a name="genomespec"> Genome specific requirements

**Download the FASTA file of your genome:** The first for every genome is a fasta file with its sequence. These can be downloaded through the 
[UCSC downloads](http://hgdownload.cse.ucsc.edu/downloads.html) page. Just select the genome you are using then click on "Full data set". Now scroll to 
the bottom of the page and there should be a file named *genomename*.fa.gz. For example, the Dolphin genome's fasta file is named *turTru2.fa.gz*. 
Sometimes your genome will be in *.2bit* format rather than *.fa*. If this is the case, see [Instruction for converting twoBitToFa](https://genome.ucsc.edu/goldenpath/help/twoBit.html)

Since there are a wide array of genomes, not every species will be included at the UCSC site. If your genome is not on the site, you will need 
to download it a different way. Unfortunately, that is outside the scope of these instructions.

**Create a BS-Seeker2 index from the fasta file:** Instructions for this can be found at [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) 
under the *bs_seeker2-build.py* section. Once you have installed everything and put BS_Seeker2 in path, the command should look something like this:

```
bs_seeker2-build.py genome.fa
```

There are a lot of defaults that can be changed so please read the [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) instructions for more details.

## <a name="Commands"> Commands </a>

Each command is accessed by typing *wgbs_tools* followed by the command name. Syntax and brief explanations are provided
using the *--help* option. For instance:

   ```
   wgbs_tools roi --help
   ```

This manual is meant to provide additional information and context for the commands. It is assumed that the user typed
*--help* for a given command before reading this manual. If there is still confusion, please email the author
(contact information at the top).

### <a name="lcommands"> Lite Commands </a>

This group of commands use the Percent Methylation BED files to create useful tables or convert them into other widely used formats.

#### <a name="bam2pm"> bam2pm </a>

Converts a bam file that was created from BS-Seeker2 to a folder of percent methylation bed files (one for each chromosome).

#### <a name="roi"> roi </a>

Finds methylation over Regions of Interest (ROIs). The inputs:

1. INPUT_TSV: A tab separated file indicating sample names and locations. Each line represents a different sample and has two columns:
  1. The first tab should be sample name.
  1. The second tab should be the path to the prefix of bed file.
1. OUT_TABLE: Name of the output table which contains the percent methylation information.
1. ROI_FILE: GTF or BED file indicating the ROI (Regions of Interest).

This command is useful if you have a few regions that you want to get the average methylation over. Some examples: gene body, 
promoters, CpG islands. A few notes:

- Since some statistical analyses would rather have counts of methylated and total reads, that information can be provided using 
the *--raw-data* option.
- You may want to mask out some regions in your analysis using the *--mask* option. For instance, CpG islands have a high concentration 
of hypomethylated CpGs, which could skew your results.

#### <a name="window"> window </a>

This command finds methylation over windows and is very similar to [roi](#roi). The main difference is that the user provides parameters
to create the regions of interest, in the form of non overlapping windows. It is highly suggested mask out CpG Islands when using this 
command due to their hypomethylation vastly skewing windows.

#### <a name="pm_stats"> pm_stats </a>

Determines basic stats of multiple pm_bed files. This will print a table where each row represents a different pm_bed file. There are 5 columns:

1. name: Unique name of file (string between prefix and suffix)
1. percentage: Percentage methylation of file
1. methylated_reads: Amount of methylated reads in pm file
1. total_reads: Amount of total reads in pm file
1. cpg_count: Number of CpGs the pm file has information for

There are a few uses for determining these basic stats:

1. Determine average methylation across a chromosome.
1. Estimate conversion efficiency. This can be done by running it on the conv_eff chromosome (see below for explanation). Then 
subtract the percentage from 1 (ie: 1 - percentage). There are 3 possible conv_eff chromosomes:
  1. If you spiked your samples with lamda DNA (known to be 100% unmethylated) prior to bisulfite conversion, run this on that 
  chromosome. Note, you will need to ensure you added the lamda sequence to info.yaml and the BS_Seeker index prior to alignment.
  1. chrM. This is the most likely choice if you didn't spike your samples.
  1. CH pm_bed files. Create pm_bed files looking at CH methylation and then run this command on all of those files. This assumes 
  there is no CH methylation (which is known to exist in plants and certain mammalian tissues).


#### <a name="pm2bg"> pm2bg </a>

Converts multiple percent methylation bed files to a single bedGraph file. This is useful for creating a trackhub. 

Track hubs, while useful for visualizing your data, require many steps to create. Briefly:

1. Convert pm_bed files to bedGraph using *wgbs_tools pm2bg*.
1. Convert the resulting bedGraph files to BigWig format using *bedGraphToBigWig*. You may need to download the converter from [UCSC](http://hgdownload.soe.ucsc.edu/admin/exe/). 
The Linux version can be found at [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig).
1. You will then need to add this information to the hub. See [Basic Track Hub Quick Start Guide](https://genome.ucsc.edu/goldenpath/help/hubQuickStart.html) 
for more information.

#### <a name="pm2dss"> pm2dss </a>

Converts multiple percent methylation bed files to DSS format. These file are necessary as input for many R packages dealing 
with methylation (such as the R packages *bsseq* and *DSS*).


#### <a name="add_genome"> add_genome </a>

Adds genome information to info.yaml file. By default it appends the *info.yaml* file in the main directory. However, you can change 
this using the *--infoyaml* option. See [TUTORIAL/README.md](https://github.com/kwdunaway/WGBS_Tools/blob/master/TUTORIAL/README.md) 
for a tutorial and examples, or see [info.yaml](#infoyaml) for more information about the *info.yaml* file.

### <a name="fcommands"> Full Commands </a>

This group of commands are all related to processing a FASTQ (or pair of FASTQ) file(s) into creating a set of Percent Methylation BED files.

#### <a name="process_se"> process_se </a>

The main pipeline to process a single end WGBS experiment. It takes in a single fastq file and outputs the following:

1. Bam containing all of the aligned reads in BS_Seeker format. See [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) for more details.
1. Folder containing a single Percent Methylation bed file for each chromosome. See [pm_bed](#pm_bed) for more information about the format.
1. Log of BS_Seeker2 run which has useful overall information about the sequencing run.
1. Bisulfite conversion efficency information.
1. Summary file containing stats about alignment and conversion efficency.

The defaults are set for determining CpG methylation in a hg38 experiment which used 100bp Illumina reads sequenced on the HiSeq 2000. 
However, the options can modified to work for your specific genome/experiment.

#### <a name="process_pe"> process_pe </a>

The main pipeline to process a paired end WGBS experiment. It takes in a pair of fastq files and outputs the same as [process_se](#process_se).

#### <a name="trim_pefq"> trim_sefq </a>

Prepares a single end fastq file for alignment by:

1. Filtering out any reads that do not pass Illumina's quality check (if that information is available in the header line of the read).
1. Trims off adapter sequence. This is done by searching the read for the adapter sequence (can be adjusted using the *--adapter* option), and 
removing all bases starting at that sequence and on.
1. Removes 10bp (can be adjusted using teh *--chew* option)from the 3' end of all reads (after adapter trimming). 

#### <a name="trim_pefq"> trim_pefq </a>

Does the same thing as [trim_sefq](#trim_sefq) except it works on a pair of fastq file. It is assumed that the reads for both files are in 
paired order (ie: the first read in the F fastq corresponds to the first read in the R fastq)

#### <a name="sumlogs"> sumlogs </a>

Summarizes BS-Seeker2 logs. This is useful to parse out the most useful information within the very long BS-Seeker2 log file.

## <a name="FileFormats"> File Formats </a>

### <a name="infoyaml"> info.yaml </a>

This is a yaml file that contains system and genome specific information. It is outside of the 
code so that any user can easily modify it to work for whatever WGBS experiment they desire. You can 
create experiment specific *.yaml* files and point to each differently using the *--infoyaml* option.
Alternatively, you can have a single file with multiple genomes, depending on your preferences.

The first two lines of *info.yaml* ***will most likely not need to be modified***:

**adapter:** The adapter sequence for single end reads. This is used to trim reads that have 
adapter contamination. It may need to be changed if you are not using Illumina sequencing.

**bs2_path:** The path for BS-Seeker2 aligner. It is set assuming BS-Seeker2 is in PATH. 
You can change this to point to a specific instance of BS Seeker 2 or if you did not put it in PATH.

The rest of the file contains genome specific information. You can have multiple genomes represented 
in the same file, but each genome must have a unique name. In order to add your genome of interest, 
use the [add_genome](#add_genome) command. The format is:

**genome:** Genome name. This needs to be unique throughout the *info.yaml* file.

**index:** Path to the BS-Seeker2 index for this genome. (2 spaces preceding) 

**fasta:** Location of a single fasta file containing all chromosomal sequences. (2 spaces preceding)

**chroms:** (2 spaces preceding, nothing after the **:** on this line)

***chromname: size*** Chromosome name (as defined by fasta/bs2index) and then size of chromosome 
in bp. Each line represents a different chromosome. You don't need every chromosome listed here, 
just the ones you want analyzed. For instance, there could be unlocalized contig names that you may 
not want to keep in your analysis. (6 spaces preceding)

You can modify the provided *info.yaml* file, or create a new one. Any commands that use this can 
point to a different file using the *--infoyaml* option.

### <a name="pm_bed"> pm_bed </a>

Percent Methylation bed (pm_bed) files contain base specific methylation information (usually for 
a single chromosome). It is based on the 9 column bed format and can be loaded into most genome browsers 
(like the [UCSC browser](http://genome.ucsc.edu/)):

1. chromosome
1. start of C
1. end of C (1 base away)
1. percent methylation and number of reads contributing to the permeth call (separate by a *-*)
1. 0 (placeholder for bed format)
1. strand (if CpG then it says + but is really a combination of both strands)
1. 0 (placeholder for bed format)
1. 0 (placeholder for bed format)
1. color (in RRR,GGG,BBB format)

***Note:*** The first line of the file is usually a header with meta information.

### <a name="meth_table"> meth_table </a>

Methylation tables are created as a way to determine the average methylation over a given area. They
are produced by both the [roi](#roi) and [window](#window) commands. The format for these files are:

1. chromosome
1. start
1. end
1. sample 1 methylation %
1. sample 2 methylation % ... (etc.)

There is also the option for printing the "raw" information. This is useful if you have further
calculations that you want to do and need read counts rather than a percentage of methylation
averaged across the region of interest. The format for these files are:

1. chromosome
1. start
1. end
1. sample 1 methylated read count
1. sample 1 total read count
1. sample 2 methylated read count
1. sample 2 total read count ... (etc.)


## <a name="statanalysis"> Statistical Analyses </a>

WGBS_Tools is meant to be a utility to take raw WGBS data and output it into a user friendly format. You can load the `pm_bed` data on a 
[genome browser](genome.ucsc.edu), get summary data, or even load the tables into R/Excel.

Unfortunately, due to its experiment specific nature, statistical analyses is outside the scope of this toolkit. However, we can give a few suggestions:

- Correct for multiple hypotheses. This is usually done with [Bonferroni](https://en.wikipedia.org/wiki/Bonferroni_correction) or 
[Benjamini](https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg_procedure) correction.
- When windowing, consider a clustering algorithm. Read [Dunaway (Cell Reports 2016)](http://www.sciencedirect.com/science/article/pii/S221112471631631X) 
for an example.
- If you want to find small (< 1000bp) DMRs, consider using the R packages DSS or bsseq. You can convert using `pm2dss`.


## <a name="VersionHistory"> Version History </a>
__0.1__:
Initial conversion of scripts to python. Improvements in multiprocessing, installation, and readability were implimented as well.

__0.0.perl__:
Most of these scripts were originally written in Perl and used in the publications:
http://www.cell.com/cell-reports/fulltext/S2211-1247(16)31631-X
https://www.ncbi.nlm.nih.gov/pubmed/28032673
However, many improvements in performance, installation ease, and understandability were implemented since this version. If you still 
want those scripts, please go checkout the branch: `git checkout perl_code`.
