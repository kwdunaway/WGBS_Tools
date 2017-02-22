# WGBS_Tools

WGBS_Tools is a versatile toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data. It is optimized for easy installation rather than efficiency. Any questions may be directed to Keith Dunaway (kwdunaway@ucdavis.edu).

### Table of Contents

1. [Installation](#Installation)
1. [Commands](#Commands)
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
**WGBS_Tools/info.yaml** = change line after bs2_path: with BS-Seeker2 path
1. Test to ensure everything installed correctly:
```
pytest tests
```



## <a name="VersionHistory"> Version History </a>
__0.1__:
Initial conversion of scripts to python. Improvements in multiprocessing, installation, and readability were implimented as well.

__0.0.perl__:
Most of these scripts were originally written in Perl and used in the publications:
http://www.cell.com/cell-reports/fulltext/S2211-1247(16)31631-X
https://www.ncbi.nlm.nih.gov/pubmed/28032673
However, many improvements in performance, installation ease, and understandability were implimented since this version. If you still want those scripts, please go to folder: WGBS_Tools/perl

