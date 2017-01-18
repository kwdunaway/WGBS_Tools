"""
Click wrapper to run all of the wgbs_tools much easier and user friendly.
"""

import os
import subprocess
import tempfile
import shutil
import click
import multiprocessing
from multiprocessing import Process
import yaml
import logging
logger = logging.getLogger(__name__)

from wgbs_tools import fastqtools
from wgbs_tools import bsseeker

NUM_CPUS = multiprocessing.cpu_count()

@click.group()
def cli():
    """A versatile toolkit to manipulate and analyze WGBS data"""
    pass


@cli.command()
@click.option('--out_dir', type=click.STRING,
              default='',
              help='Directory to put all outfiles. Default <EMPTY> (Creates '
                   'directory with same name as out_prefix)')
@click.option('--genome', type=click.STRING,
              default='hg38',
              help='Genome used for alignment and analysis. Default: hg38')
@click.option('--noadap-bs2params', 'noadap_bs2_params', type=click.STRING,
              default='-m 3 -f bam ',
              help='Parameters passed to BS Seeker 2 for alignment of reads '
                   'without adapter contamination. Default: -m 3 -f bam')
@click.option('--adaptrim-bs2params', 'adaptrim_bs2_params', type=click.STRING,
              default='-m 2 -f bam ',
              help='Parameters passed to BS Seeker 2 for alignment of reads '
                   'with adapter contamination trimmed out. Default: -m 2 -f '
                   'bam')
@click.option('--trimmed/--not-trimmed',
              default=False,
              help='Input fastq file is already trimmed for adapter sequence.'
                   'Default: --not-trimmed')
@click.option('--cgibed/--no-cgibed',
              default=True,
              help='For each PerMeth bed file, creates a separate one without '
                   'CPG island information in it. This is useful for looking '
                   'at windows without CGI bias. Default: --cgibed (creates '
                   'the separate permeth bed files)')
@click.option('--threads', type=click.INT,
              default=NUM_CPUS,
              help='Number of threads used when multiprocessing. '
                   'Default: Number of system CPUs')
@click.option('--working_dir', type=click.STRING,
              default='',
              help='Working directory where temp files are written and read. '
                   'Default: <EMPTY> (Uses tempfile.mkdtemp() to define a '
                   'temperary working directory)')
@click.option('--infoyaml', type=click.STRING,
              default='',
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. Default: info.yaml')
@click.argument('input', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def main_pipeline(in_fastq, out_prefix, out_dir, genome, noadap_bs2_params,
                  trimmed, cgibed, threads, working_dir, infoyaml):
    """
    Main fastq to permeth bed pipeline.

    #!/bin/bash
#
#SBATCH --workdir /share/lasallelab/CM_WGBS_CordBlood_Batch2/
#SBATCH -c 12                               # number of processors
#SBATCH -N 1                                # number of nodes

PATH=$PATH:/share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/
module load bowtie/1.1.1
module load samtools/0.1.19
module load sratoolkit/2.4.2-3
module load bedtools2/2.25.0
export PYTHONPATH=/share/lasallelab/pysam/lib/python2.7/site-packages/

#Complete Run in hg38 for JLCM007A
gunzip -c raw_sequences/JLCM007A*fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > raw_sequences/JLCM007A_filtered.fq
gzip raw_sequences/JLCM007A_filtered.fq
perl /share/lasallelab/programs/perl_script/adapter_split.pl raw_sequences/JLCM007A_filtered.fq.gz raw_sequences/JLCM007A_noadap.fq.gz raw_sequences/JLCM007A_withadap.fq.gz
perl /share/lasallelab/programs/perl_script/adapter_trimmer.pl raw_sequences/JLCM007A_withadap.fq.gz raw_sequences/JLCM007A_trimmed.fq.gz 45 10
mkdir JLCM007A
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3 -f bam -g /share/lasallelab/genomes/hg38/hg38.fa -d /share/lasallelab/genomes/hg38/BSseek2_refgen/ -i raw_sequences/JLCM007A_noadap.fq.gz -o JLCM007A/JLCM007A_noadap.bam
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 80 -m 2 -f bam -g /share/lasallelab/genomes/hg38/hg38.fa -d /share/lasallelab/genomes/hg38/BSseek2_refgen/ -i raw_sequences/JLCM007A_trimmed.fq.gz -o JLCM007A/JLCM007A_trimmed.bam
samtools sort JLCM007A/JLCM007A_noadap.bam JLCM007A/JLCM007A_noadap_sorted
samtools sort JLCM007A/JLCM007A_trimmed.bam JLCM007A/JLCM007A_trimmed_sorted
samtools merge JLCM007A/JLCM007A.bam JLCM007A/JLCM007A_noadap_sorted.bam JLCM007A/JLCM007A_trimmed_sorted.bam
samtools view JLCM007A/JLCM007A.bam > JLCM007A/JLCM007A.sam
mkdir JLCM007A/tmp
perl /share/lasallelab/programs/perl_script/SAMsorted_to_permeth.pl JLCM007A/JLCM007A.sam JLCM007A/tmp/PerMeth_JLCM007A JLCM007A hg38 CG combined 1
mkdir JLCM007A/PerMeth_JLCM007A
# Unnecessary because incorporated into previous script
#perl /share/lasallelab/programs/perl_script/gbcompliance.pl hg38 JLCM007A/tmp/PerMeth_JLCM007A_ JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_ JLCM007A JLCM007A
rm -r JLCM007A/tmp
mkdir JLCM007A/NoCGI_Permeth_JLCM007A
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr1.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr1.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr2.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr2.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr3.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr3.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr4.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr4.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr5.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr5.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr6.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr6.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr7.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr7.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr8.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr8.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr9.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr9.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr10.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr10.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr11.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr11.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr12.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr12.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr13.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr13.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr14.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr14.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr15.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr15.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr16.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr16.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr17.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr17.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr18.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr18.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr19.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr19.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr20.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr20.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr21.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr21.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chr22.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chr22.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chrX.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chrX.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chrY.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chrY.bed
bedtools subtract -a JLCM007A/PerMeth_JLCM007A/PerMeth_JLCM007A_chrM.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > JLCM007A/NoCGI_Permeth_JLCM007A/NoCGI_Permeth_JLCM007A_chrM.bed

    """
    if working_dir == '':
        workingdir = tempfile.mkdtemp()
    else:
        workingdir = working_dir
    #Load data
    stream = file(infoyaml, 'r')
    info_dict = yaml.safe_load(stream)
    bs2_path = info_dict['bs2_path']
    bs2_index = info_dict[genome]['bs2_index']
    fasta = info_dict[genome]['fasta']
    chroms = info_dict[genome]['chroms']

    #Name temp files
    qualfil_fastq = os.path.join(workingdir, '_filtered.fq'.format(out_prefix))
    noadap_fq = os.path.join(workingdir, '_noadap.fq'.format(out_prefix))
    adaptrim_fq = os.path.join(workingdir, '_trimmed.fq'.format(out_prefix))
    noadap_bam = os.path.join(workingdir, '_noadap.bam'.format(out_prefix))
    adaptrim_bam = os.path.join(workingdir, '_adaptrim.bam'.format(out_prefix))
    noadap_sorted = os.path.join(workingdir, '_noadap_sorted'
                                 .format(out_prefix))
    adaptrim_sorted = os.path.join(workingdir, '_adaptrim_sorted'
                                   .format(out_prefix))
    full_bam = os.path.join(workingdir, '.bam'.format(out_prefix))
    full_sam = os.path.join(workingdir, '.sam'.format(out_prefix))

    #Create output directory
    if out_dir == '':
        out_dir = out_prefix
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #Filter fastq file
    logging.info('Filtering out quality failed reads from fastq file')
    fastqtools.qual_filter_fastq(in_fastq, qualfil_fastq)

    #Remove adapter contamination from reads
    logging.info('Removing adapter contamination and split fastq files')
    fastqtools.adapter_remove(qualfil_fastq, noadap_fq, adaptrim_fq,
                             info_dict['adapter'])

    #Align to genome
    logging.info('Aligning noadap contamination reads to {}'.format(genome))
    noadap_bs2_params = '--bt-p {} {}'.format(threads, noadap_bs2_params)
    bsseeker.align_bs2(bs2_path, noadap_bs2_params, fasta, bs2_index, noadap_fq,
                       noadap_bam)
    logging.info('Aligning adaptrim contamination reads to {}'.format(genome))
    adaptrim_bs2_params = '--bt-p {} {}'.format(threads, adaptrim_bs2_params)
    bsseeker.align_bs2(bs2_path, adaptrim_bs2_params, fasta, bs2_index,
                       adaptrim_fq, adaptrim_bam)
    #TODO: Combine alignment logs

    #Sort bam files
    command = 'samtools sort -@ {} {} {}'\
        .format(threads, noadap_bam, noadap_sorted)
    logging.info(command)
    subprocess.check_call(command, shell=True)
    command = 'samtools sort -@ {} {} {}'\
        .format(threads, adaptrim_bam, adaptrim_sorted)
    logging.info(command)
    subprocess.check_call(command, shell=True)

    #Merge bam files
    command = 'samtools merge -@ {} {} {}.bam {}.bam'\
        .format(threads, full_bam, noadap_sorted, adaptrim_sorted)
    logging.info(command)
    subprocess.check_call(command, shell=True)

    #Convert to sam format
    command = 'samtools view {} > {}'.format(full_bam, full_sam)
    logging.info(command)
    subprocess.check_call(command, shell=True)

    #Convert sam to permeth bed files (percent methylation bed files)
    #This also only prints CpGs on chromosomes in ranges defined in the yaml

    #Creates no CGI permeth bed files (if applicable)
#    if cgibed:






