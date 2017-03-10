# WGBS_Tools Examples

This will go through examples of how to use wgbs_tools. All of the data is provided, however it is necessary you follow the installation instructions before proceeding.

### Table of Contents

1. [Workflow](#Workflow)
  1. [Prerequisites](#Prerequisites)
  1. [Instructions](#Instructions)
1. [Commands](#Commands)
  1. [Adjust existing files](#acommands)

## <a name="Workflow"> Workflow </a>

Here is a basic workflow (the wgbs_tools commands are in red):

![work_outline](work_outline.png)


## Load genomic information

Before processing your sample data, you first need to ensure you have the genomic specific information accessible. That means:

1. Download the FASTA file of your genome. This can be done through the [UCSC downloads](http://hgdownload.cse.ucsc.edu/downloads.html) page. Just select
the genome you are using then click on "Full data set". Now scroll to the bottom of the page and there should be a file named *genomename*.fa.gz. For example,
the Dolphin genome's fasta file is named *turTru2.fa.gz*.

1. Create a BS-Seeker2 index from the fasta file. Instructions for this can be found at [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) 
under the *bs_seeker2-build.py* section.

You then need to add the genome specific information to you *info.yaml* file. You can do this one of two ways:

1. Add your genome to the *info.yaml* file by using the following command:

```
wgbs_tools 
```

1. 