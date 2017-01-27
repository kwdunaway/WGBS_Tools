# WGBS_Tools

WGBS_Tools is a versatile toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data. It is optimized for easy installation rather than efficiency. Any questions may be directed to Keith Dunaway (kwdunaway@ucdavis.edu).

### Table of Contents

1. [Installation](#Installation)

1. [Pipelines](#Pipelines)

1. [Commands](#Commands)

1. [Version History](#VersionHistory)

## <a name="Installation"> Installation </a>

The following programs need to be installed before you can use WGBS_Tools:

  [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml)
  [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)
  [Pysam](https://github.com/pysam-developers/pysam)
  [Bedtools2](https://github.com/arq5x/bedtools2)
  [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/)

### python packages:
pyyaml

### path changes:
Since every system has a different path, you need to edit the paths in the following files:

1. **WGBS_Tools/info.yaml** = change line after bs2_path: with BS-Seeker2 path


## <a name="VersionHistory"> Version History </a>
__0.1__:
Initial conversion of scripts to python. Improvements in multiprocessing, installation, and readability were implimented as well.

__0.0.perl__:
Most of these scripts were originally written in Perl and used in the publications:
http://www.cell.com/cell-reports/fulltext/S2211-1247(16)31631-X
https://www.ncbi.nlm.nih.gov/pubmed/28032673
However, many improvements in performance, installation ease, and understandability were implimented since this version. If you still want those scripts, please go to folder: WGBS_Tools/perl

