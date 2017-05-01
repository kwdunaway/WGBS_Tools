# WGBS_Tools TUTORIAL

This will go through examples of how to use wgbs_tools. All of the data is provided, however 
it is necessary you follow the installation instructions before proceeding.

**READ BEFORE STARTING TUTORIAL:** Everything is assumed you are running these commands in the 
TUTORIAL folder. So change your directory to this folder:

```
cd TUTORIAL
```

### Table of Contents

1. [Lite Workflow](#lWorkflow)
    1. [Lite version: Add a genome to .yaml ](#lAddGenome)
    1. [Convert bam to pm_bed](#bam2pm)
    1. [Analyze your experiment](#analysis)
        1. [Windowing](#window)
        1. [Regions of Interest](#roi)
    1. [Convert pm_bed to another format](#convert)
        1. [pm_bed to DSS](#expm2dss)
1. [Full Workflow](#fWorkflow)
    1. [Full version: Add a genome to .yaml ](#fAddGenome)
    1. [Process fastq file in single step](#processse)
    1. [Process fastq file in multiple steps](#processstep)
        1. [Filter and trim reads](#adaptrim)
        1. [Align using BS Seeker2](#alignse)
        1. [Convert bam to pm_bed](#fbam2pm)


## <a name="Workflow"> Lite Workflow </a>

Here is a basic workflow of the Lite commands (the wgbs_tools commands are in red):

![work_outline](lite_outline.png)

### <a name="lAddGenome"> Lite version: Add a genome to .yaml </a>

The first thing you need to do before starting a new experiment is to add the reference genome information to the .yaml file
if it is not already there. For real genomes, you will need to download the fasta file and create a BS-Seeker2 index.
Instructions for that can be found at [Full version: Add a genome to .yaml](#fAddGenome) but for right now we will
only focus on this tutorial genome. Also, you will usually want to edit the *info.yaml* file in the root directory of WGBS_Tools
but for this tutorial we will be editing and using the *example.yaml* file.

The .yaml file contains genome specific information necessary for most processes in WGBS_Tools. It is designed for a user
to easily add more genomes through the add_genome command. First though, take a look at the format of the file:

```
less example.yaml
```

The first two lines of the file contain information necessary for functions in the full version of WGBS_Tools, so we will
ignore those for now. The next line designates a genome name. Everything that is indented after that is information pertaining
to that genome. For more information on the *.yaml* format, see the https://github.com/kwdunaway/WGBS_Tools#infoyaml page.

Now, try adding hg38 genomic information to the file then look at the difference it makes:

```
wgbs_tools add_genome --infoyaml example.yaml --force hg38 /path/to/fasta.fa /path/to/indexdir/
less example.yaml
```

You will notice that chromosome information about hg38 was added to the file (the expected 1-22, X, Y, M). However, these are only the main
chromosomes. If you want the extra chromosomes that include "_" in the name (for example: *_alt* and *_random* chromosomes) use the *--all* option.

You can check to see if your command worked correctly by looking at the expected file(s) in *correct_output*. This is true for all examples, 
in this tutorial. You can manually inspect the differences to see if everything is correct, or **diff** the files. For example:

```
diff example.yaml correct_output/example.yaml
```

### <a name="bam2pm"> Convert bam to pm_bed </a>

Both BS Seeker2 and Bismark output their data to sam/bam format by default. The methylation call strings are slightly 
different between them but `bam2pm` can distinguish these differences and works for both formats. To keep things simple,
everything in this tutorial will be derived from BS-Seeker2 formats. We can convert BS-Seeker2 bam file to pm_bed format by:

```
wgbs_tools bam2pm --genome example_genome --infoyaml example.yaml ex.bam out
```

### <a name="analysis"> Analyze your experiment </a>

You may have a variety of different experiments, but they usually all boil down to the question: "Are there areas of the 
genome that are hypo/hyper methylated?". To answer this, you need to compare your samples against each other (specifically, 
control samples to the experimental ones). But first you need to create a file that identifies where all of your samples 
are and what sample they came from. Look at `example_input.tsv` in this folder for an example of what this tab separated file 
(tsv) looks like. It has two columns:

1. Sample name. This can be anything but I suggest you keep it short for easier readability
1. Prefix to pm_bed files. When creating pm_bed files, they are broken up by chromosome. This is the path all the way to that prefix.

You will notice that the path for all of the samples leads to the `data/pmbeds` folder. If you look in the folder, each sample 
only has one `pm_bed` file (for chrT). The chromosomes are defined by the `.yaml` file and `example.yaml`, which has only one 
chromosome (chrT) for the genome `example_genome`. In a real genome there will be multiple pm_bed files corresponding to the 
multiple chromosomes of the genome.

You will also notice that there is one zipped file in the `pmbeds` folder. This is to demonstrate that all of these commands 
work on zipped or unzipped pm_bed files.

#### <a name="window"> Windowing </a>

**Note:** Windowing the genome takes a while to run due to the large size of most genomes. Expect it to take hours.

The `window` command breaks the genome up into equal, non-overlapping sections (aka windows). The default size is 20000bp with 
at least 20 CpGs represented in every sample, which can all be changed with command options. An example of how to use the command:

```
wgbs_tools window --window-size 10 --min-cpgs 1 --infoyaml example.yaml --genome example_genome data/example_input.tsv windowed.txt
```

If you look at the `example.yaml`, chrT is set at 100bp. So, you would expect there to be windows 10 windows at 10bp each. However, 
you will notice that in the windowed.txt file, there are only 4 windows. This is because if there is no data on a window for any of 
the samples, that window is not printed.

To get a more complete idea of what's going on, let's look at the raw data:

```
wgbs_tools window --raw-data windowed_rawdata.txt --window-size 10 --min-cpgs 1 --infoyaml example.yaml --genome example_genome data/example_input.tsv windowed.txt
```

This is the same command and the previous, with just one additional option: `--raw-data windowed_rawdata.txt`. If you look at the 
file created, you can see the number of methylated and total reads are reported, instead of just the average methylation across 
the window. Some statistical analyses would rather have this raw data as input.

Now try the command:
 
```
wgbs_tools window --window-size 10 --min-cpgs 2 --infoyaml example.yaml --genome example_genome data/example_input.tsv windowed2.txt
```

If you compare the `windowed.txt` and `windowed2.txt`, you will notice that there are 4 lines in `windowed.txt` and only 1 in 
`windowed2.txt`. By changing the minimum number of CpGs required per window, it will remove any windows that do not meet this requirement. 

It is highly suggested that when you window the genome, you mask out the CpG Islands. This is because they are both highly concentrated 
in CpG and unusually hypomethylated compared to the rest of the genome. To do this, download the CpG island locations and then use the 
`--mask` option:

```
wgbs_tools window --mask data/mask.bed --window-size 10 --min-cpgs 1 --infoyaml example.yaml --genome example_genome data/example_input.tsv windowed3.txt
```

The `--mask` option is also useful if you want to mask out other areas of the genome.

#### <a name="roi"> Regions of Interest </a>

There will be areas of the genome that you are interested in analyzing. For instance, you may want to see if there is hypomethylation 
in the gene body. In that case, you will need to get a bed or gtf file of these areas. If you want to download annotation files for your 
genome, try [UCSC downloads](http://hgdownload.cse.ucsc.edu/downloads.html) or [UCSC tables](http://genome.ucsc.edu/cgi-bin/hgTables?command=start).

This example will take the `roi.bed` as a region of interest file.

```
wgbs_tools roi --raw-data roi_table_rawdata.txt --min-cpgs 1 data/example_input.tsv roi_table.txt data/roi.bed
```

Many of options in `window` are also available in `roi`:

```
wgbs_tools roi --raw-data roi_table_rawdata.txt --min-cpgs 1 data/example_input.tsv roi_table.txt data/roi.bed
wgbs_tools roi --min-cpgs 2 data/example_input.tsv roi_table2.txt data/roi.bed
wgbs_tools roi --mask data/mask.bed --min-cpgs 1 data/example_input.tsv roi_table3.txt data/roi.bed
```

One important difference is that `roi` does not need a `.yaml` file to run. This is because while `window` needs chromosome lengths to 
subset the genome, roi only uses the areas designated in the `roi.bed` file.

### <a name="convert"> Convert pm_bed to another format </a>

You may want to convert your data to another format. Some of the most common formats have converters built into `wgbs_tools`.

#### <a name="expm2dss"> pm_bed to DSS </a>

If you want to find small (< 1000bp) DMRs you should consider using the R packages DSS or bsseq. In order to use these modules, you 
need to convert your data to DSS format:

```
mkdir DSS
wgbs_tools pm2dss data/pmbeds/pm DSS/pm
```

## <a name="Workflow"> Full Workflow </a>

Here is a basic workflow of the Full commands (the wgbs_tools commands are in red):

![work_outline](work_outline.png)

## <a name="fAddGenome"> Full version: Add a genome to .yaml </a>

Before processing your sample data, you first need to ensure you have the genomic specific information accessible. That means:

1. Download the FASTA file of your genome. This can be done through the [UCSC downloads](http://hgdownload.cse.ucsc.edu/downloads.html) page. 
Just select the genome you are using then click on "Full data set". Now scroll to the bottom of the page and there should be a file named 
*genomename*.fa.gz. For example, the Dolphin genome's fasta file is named *turTru2.fa.gz*.

1. Create a BS-Seeker2 index from the fasta file. Instructions for this can be found at [BS-Seeker2](https://github.com/BSSeeker/BSseeker2) 
under the *bs_seeker2-build.py* section.

You then need to add the genome specific information to you *info.yaml* file. You can do this one of two ways:

1. Add your genome to the *info.yaml* file by using the *add_genome* command. Look at the example.yaml file before and after running this 
command to see how *add_genome* changes the file.

```
less example.yaml
wgbs_tools add_genome --infoyaml example.yaml --force hg38 /path/to/fasta.fa /path/to/indexdir/
less example.yaml
```

You will notice that only the main chromosomes were added to the file (the expected 1-22, X, Y, M). If you want the extra chromosomes that 
include "_" in the name (for example: *_alt* and *_random* chromosomes) use the *--all* option.

For this and all other examples, you can check to see if your command worked correctly by looking at the corresponding file(s) in 
*correct_output*. You can manually inspect them or **diff** the files. For example:

```
diff example.yaml correct_output/example.yaml
```

## <a name="processse"> Process fastq file in single step </a>

Once the genome/system specific information is in the *.yaml* file, you are ready to start processing your experiment. 
An example of how to process a single-end WGBS experiment in a single command:

```
wgbs_tools process_se --genome example_genome --infoyaml example.yaml data/example.fq out
```

This should produce the directory `out`. It should look just like the one in the `correct` folder. This is the end product of `process_se`. It contains:

- bam file containing all alignment results (both from trimmed and noadap)
- bam index of the bam file (.bai)
- summary.txt file containing a summary of some useful statistics of the alignment
- permethbed folder containing a bed file for each chromosome. Since this test set only has chrT, it only created 1 bed file.

If you want to see the intermediate files, try running the command:

```
rm -r out
wgbs_tools process_se --genome example_genome --infoyaml example.yaml --working-dir workingdir data/example.fq out
```

This will produce an extra folder called `workingdir` and contains all of the intermediate files. The default settings remove these files after processing but if you want to keep them then you should use the `--working-dir` option. Once you are done looking at these file, removed them:

```
rm -r out
rm -r workingdir
```
## <a name="processstep"> Process fastq file in multiple steps </a>

While `process_se` is useful if you want to use standard settings, you could also process the fastq to pm_bed in individual steps. This is useful if you want to customize certain parameters that are not available in the `process_se` command.

### <a name="adaptrim"> Filter and trim reads </a>

The first step is to quality filter and adapter trim the reads:

```
wgbs_tools trim_sefq data/example.fq out
```

This should produce 3 files:

1. `out_filtered.fq.gz` which is a file containing reads that passed Illumina's quality filter.
1. `out_noadap.fq.gz` containing Illumina's quality filtered reads which did not have adapter contamination. They are chewed back 10bp.
1. `out_trimmed.fq.gz` containing Illumina's quality filtered reads which had adapter contamination. The adapter sequece was removed 
then the reads were chewed back 10bp.

It has been found that the bases closest to the 3' adapter had consistantly lower quality and gives incorrect methylation calls. So, 
the `chew` was implimented to remove those low quality bases. The `chew` length can be changed although 10bp is consistant with the literature.

Also, a minimum read length of 35bp is default but that can be changed with the `min-read` option.

### <a name="alignse"> Align using BS Seeker2 </a>

There are currently two main aligners for WGBS data: BS Seeker2 and Bismark. This example will only use BS Seeker2. Since these examples 
are going step by step, the alignment step is done by directly calling the program.

However, before alignment, you will need to create a genome index from the fasta file:

```
bs_seeker2-build.py -f example_genome/genome.fa -d new_ref_genome
```

Now you can align your fastq files to the genome:

```
bs_seeker2-align.py --bt-p 8 -m 3 -f bam -g example_genome/genome.fa -d new_ref_genome/ -i out_noadap.fq.gz -o out_noadap.bam
bs_seeker2-align.py --bt-p 8 -m 2 -f bam -g example_genome/genome.fa -d new_ref_genome/ -i out_trimmed.fq.gz -o out_trimmed.bam
```

By default, 2 reads are allowed as mismatches in the trimmed alignment while 3 are allowed in the noadap alignment. This is due to 
the trimmed reads being shorter than the noadap reads.

Once the alignment is complete, you will need to combine the two bam files:

```
samtools merge out.bam out_noadap.bam out_trimmed.bam
```

Finally, you need to create an index:

```
samtools index out.bam
```

### <a name="fbam2pm"> Convert bam to pm_bed </a>

At this point we have done everything before the steps already described in the Lite version. If you would like instructions on how
to proceed, please see the Lite version's [Convert bam to pm_bed](#fbam2pm) section.
