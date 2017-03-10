# WGBS_Tools Examples

This will go through examples of how to use wgbs_tools. All of the data is provided, however it is necessary you follow the installation instructions before proceeding.

Everything is assumed you are running these commands in the EXAMPLES folder. So change your directory to this folder:

```
cd EXAMPLES
```

### Table of Contents

1. [Workflow](#Workflow)
1. [Add genome to .yaml](#AddGenome)
1. [Process fastq file](#processse)

## <a name="Workflow"> Workflow </a>

Here is a basic workflow (the wgbs_tools commands are in red):

![work_outline](work_outline.png)


## <a name="AddGenome"> Add a genome to .yaml </a>

Before processing your sample data, you first need to ensure you have the genomic specific information accessible. That means:

1. Download the FASTA file of your genome. This can be done through the [UCSC downloads](http://hgdownload.cse.ucsc.edu/downloads.html) page. Just select
the genome you are using then click on "Full data set". Now scroll to the bottom of the page and there should be a file named *genomename*.fa.gz. For example,
the Dolphin genome's fasta file is named *turTru2.fa.gz*.

1. Create a BS-Seeker2 index from the fasta file. Instructions for this can be found at [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) 
under the *bs_seeker2-build.py* section.

You then need to add the genome specific information to you *info.yaml* file. You can do this one of two ways:

1. Add your genome to the *info.yaml* file by using the *add_genome* command. Look at the example.yaml file before and after running
this command to see how *add_genome* changes the file.

```
less example.yaml
wgbs_tools add_genome --infoyaml example.yaml --force hg38 /path/to/fasta.fa /path/to/indexdir/
less example.yaml
```

You will notice that only the main chromosomes were added to the file (the expected 1-22, X, Y, M). If you want the extra chromosomes that include "_" in the name
(for example: *_alt* and *_random* chromosomes) use the *--all* option.

You can check to see if this worked correctly by looking at the file that should have been created in *correct/example.yaml*. Alternatively, you can diff the files:

```
diff example.yaml correct/example.yaml
```

## <a name="processse"> Process fastq file </a>

Once the genome/system specific information is in the *.yaml* file, you are ready to start processing your experiment. An example of how to 
process a single-end WGBS experiment in a single command:

```
wgbs_tools process_se
```

