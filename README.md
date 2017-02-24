# WGBS_Tools

WGBS_Tools is a versatile toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data. It is optimized for easy installation rather than efficiency. Any questions may be directed to Keith Dunaway (kwdunaway@ucdavis.edu).

### Table of Contents

1. [Installation](#Installation)
1. [File Formats](#FileFormats)
  1. [infoyaml](#info.yaml)
  1. [pm_bed](#pm_bed)
  1. [meth_table](#meth_table)
1. [Commands](#Commands)
  1. [align](#align)
  1. [roi](#roi)
  1. [window](#window)
  1. [pm_stats](#pm_stats)
  1. [pm2bg](#pm2bg)
  1. [pm2dss](#pm2dss)
1. [Version History](#VersionHistory)

## <a name="Installation"> Installation </a>

**Note:** It is recommended that you set up a virtual environment prior
to installing wgbs_tools. The author uses *venv* and will assume the user
does the same. Using another virtual environment will most likely not make
a difference. However, other virtual environments are not supported.

### Prerequisites:

Before installation, essure that he following programs are installed:

- [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml)
- [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)
- [Pysam](https://github.com/pysam-developers/pysam)
- [Bedtools2](https://github.com/arq5x/bedtools2)
- [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/)

Also, the following python packages:

- click
- pybedtools
- pysam
- pytest
- pyyaml

### Instructions:
1. Change directory to WGBS_tools and type the following:

   ```
   python setup.py develop
   ```
1. Since every system has a different path, you need to edit the paths in the following file:

   **WGBS_Tools/info.yaml** = change line after *bs2_path:* with BS-Seeker2 path
1. Test to ensure everything installed correctly:

   ```
   pytest tests
   ```
1. Instructions for how to use the commands are below as well as accessible by typing:

   ```
   wgbs_tools --help
   ```
   
1. You can also get instructions for a specific command by typing the command name with --help. For instance:

   ```
   wgbs_tools align --help
   ```

1. The final step is to edit (or create your own) *info.yaml* file. You must add your own genomic information.
Chromosome lengths for genomes hg38, mm10, and bosTau6 are provided. However, you will still need to edit the
*index* and *fasta* locations to point to where those files are on your system. See [infoyaml](#info.yaml) for more details.


## <a name="FileFormats"> File Formats </a>

### <a name="infoyaml"> info.yaml </a>



### <a name="pm_bed"> pm_bed </a>

Percent Methylation bed (pm_bed) format is a custome format to hold base specific methylation information.
It is based on a 9 column bed format and can be loaded into most genome browsers (like the UCSC browser).
The pm_bed has 7 columns:

1. chromosome
1. start of C
1. end of C (1 base away)
1. percent methylation and number of reads contributing to the permeth call (separate by a *-*)
1. 0 (placeholder for bed format)
1. strand (if CpG then it says + but is really a combination of both strands)
1. 0 (placeholder for bed format)
1. 0 (placeholder for bed format)
1. color (in RRR,GGG,BBB format)

### <a name="meth_table"> meth_table </a>

Methylation tables are created as a way to determine the average methylation over a given area. They
are produced by both the [roi](#roi) and [window](#window) commands. The format for these files are:

bed (pm_bed) format is a custome format to hold base specific methylation information.
It is based on a 9 column bed format and can be loaded into most genome browsers (like the UCSC browser).
The pm_bed has 7 columns:

1. chromosome
1. start
1. end

## <a name="Commands"> Commands </a>

### <a name="align"> align </a>
The main alignment pipeline. It takes in a fastq file and outputs the following:
1. Bam containing all of the aligned reads in BS_Seeker format. See [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) for more details.
1. Log of BS_Seeker2 run which has useful overall information about the sequencing run.
1. You need to edit the info.yaml with your genome information (see info.yaml)


### <a name="roi"> roi </a>
  roi        Calls methylation over ROIs.

### <a name="window"> window </a>
  window     Calls methylation over windows.

### <a name="pm_stats"> pm_stats </a>
  pm_stats   Gets stats of multiple pm_bed files.

### <a name="pm2bg"> pm2bg </a>
  pm2bg      Converts pm_bed to bedgraph format.

### <a name="pm2dss"> pm2dss </a>
  pm2dss     Converts pm_bed to dss format.


## <a name="VersionHistory"> Version History </a>
__0.1__:
Initial conversion of scripts to python. Improvements in multiprocessing, installation, and readability were implimented as well.

__0.0.perl__:
Most of these scripts were originally written in Perl and used in the publications:
http://www.cell.com/cell-reports/fulltext/S2211-1247(16)31631-X
https://www.ncbi.nlm.nih.gov/pubmed/28032673
However, many improvements in performance, installation ease, and understandability were implimented since this version. If you still want those scripts, please go to folder: WGBS_Tools/perl

