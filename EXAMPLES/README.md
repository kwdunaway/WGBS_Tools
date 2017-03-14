# WGBS_Tools Examples

This will go through examples of how to use wgbs_tools. All of the data is provided, however it is necessary you follow the installation instructions before proceeding.

Everything is assumed you are running these commands in the EXAMPLES folder. So change your directory to this folder:

```
cd EXAMPLES
```

### Table of Contents

1. [Workflow](#Workflow)
1. [Add genome to .yaml](#AddGenome)
1. [Process fastq file in single step](#processse)
1. [Process fastq file in multiple steps](#processstep)
  1. [Filter and trim reads](#adaptrim)
  1. [Align using BS Seeker2](#alignse)
  1. [Convert bam to pm_bed](#bam2pm)
1. [Window the genome](#window)

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

## <a name="processse"> Process fastq file in single step </a>

Once the genome/system specific information is in the *.yaml* file, you are ready to start processing your experiment. An example of how to 
process a single-end WGBS experiment in a single command:

```
wgbs_tools process_se --genome example --infoyaml example.yaml example.fq outdir
```

This should produce the directory `outdir`. It should look just like the one in the `correct` folder. This is the end product of `process_se`. It contains:

- bam file containing all alignment results (both from trimmed and noadap)
- bam index of the bam file (.bai)
- summary.txt file containing a summary of some useful statistics of the alignment
- permethbed folder containing a bed file for each chromosome. Since this test set only has chrT, it only created 1 bed file.

If you want to see the intermediate files, try running the command:

```
rm -r outdir
wgbs_tools process_se --genome example --infoyaml example.yaml --working-dir workingdir example.fq outdir
```

This will produce an extra folder called `workingdir` and contains all of the intermediate files. The default settings remove these files after processing but
if you want to keep them then you should use the `--working-dir` option.

## <a name="processstep"> Process fastq file in multiple steps </a>

While `process_se` is useful if you want to use standard settings, you could also process the fastq to pm_bed in individual steps. This is useful
if you want to customize certain parameters that are not available in the `process_se` command.

### <a name="adaptrim"> Filter and trim reads </a>

The first step is to quality filter and adapter trim the reads:

```
wgbs_tools trim_sefq example.fq outdir
```

This should produce 3 files:

1. `outdir_filtered.fq.gz` which is a file containing reads that passed Illumina's quality filter.
1. `outdir_noadap.fq.gz` containing Illumina's quality filtered reads which did not have adapter contamination. They are chewed back 10bp.
1. `outdir_trimmed.fq.gz` containing Illumina's quality filtered reads which had adapter contamination. The adapter sequece was removed then the reads were chewed back 10bp.

It has been found that the bases closest to the 3' adapter had consistantly lower quality and gives incorrect methylation calls. So, the `chew` was implimented to
remove those low quality bases. The `chew` length can be changed although 10bp is consistant with the literature.

Also, a minimum read length of 35bp is default but that can be changed with the `min-read` option.

### <a name="alignse"> Align using BS Seeker2 </a>

There are currently two main aligners for WGBS data: BS Seeker2 and Bismark. This example will only use BS Seeker2. Since these examples are going step by step, 
the alignment step is done by directly calling the program.

However, before alignment, you will need to create a genome index from the fasta file:

```
bs_seeker2-build.py -f example_genome/genome.fa -d new_ref_genome
```

Now you can align your fastq files to the genome:

```
bs_seeker2-align.py --bt-p 8 -m 3 -f bam -g example_genome/genome.fa -d new_ref_genome/ -i outdir_noadap.fq.gz -o outdir_noadap.bam
bs_seeker2-align.py --bt-p 8 -m 2 -f bam -g example_genome/genome.fa -d new_ref_genome/ -i outdir_trimmed.fq.gz -o outdir_trimmed.bam
```

By default, 2 reads are allowed as mismatches in the trimmed alignment while 3 are allowed in the noadap alignment. This is due to the trimmed reads being
shorter than the noadap reads.

Once the alignment is complete, you will need to combine the two bam files:

```
samtools merge outdir.bam outdir_noadap.bam outdir_trimmed.bam
```

Finally, you need to create an index:

```
samtools index outdir.bam
```

### <a name="bam2pm"> Convert bam to pm_bed </a>

Both BS Seeker2 and Bismark output their data to sam/bam format by default. However, the methylation call strings are slightly different between them. 
So, use the correct converter: `bam2pm` for BS Seeker2 and `bisbam2pm` for Bismark. Working from the previous data, we can convert BS Seeker2 bam
file to pm_bed format by:

```
wgbs_tools bam2pm --genome example --infoyaml example.yaml outdir.bam outdir
```

Once in pm_bed format, you can analyze the experiment using a variety of tools in the `wgbs_tools` module.

### <a name="window"> Window the genome </a>

One question you may have is "Are there areas of the genome that are hypo/hyper methylated?". To answer this, you 