"""
Click wrapper to run all of the wgbs_tools much easier and user friendly.
"""

import os
import glob
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
from wgbs_tools import utilities
from wgbs_tools import samutils
from wgbs_tools import permethbed

NUM_CPUS = multiprocessing.cpu_count()

@click.group()
def cli():
    """A versatile toolkit to manipulate and analyze WGBS data"""
    pass


@cli.command()
@click.option('--out_dir', type=click.STRING,
              default='',
              help='Directory to put all outfiles. '
                   'Default: <current working directory>')
@click.option('--genome', type=click.STRING,
              default='hg38',
              help='Genome used for alignment and analysis. '
                   'Default: hg38')
@click.option('--noadap-bs2params', 'noadap_bs2_params', type=click.STRING,
              default='-m 3 -f bam',
              help='Parameters passed to BS Seeker 2 for alignment of reads '
                   'without adapter contamination. '
                   'Default: -m 3 -f bam')
@click.option('--adaptrim-bs2params', 'adaptrim_bs2_params', type=click.STRING,
              default='-m 2 -f bam',
              help='Parameters passed to BS Seeker 2 for alignment of reads '
                   'with adapter contamination trimmed out. '
                   'Default: -m 2 -f bam')
@click.option('--trimmed/--not-trimmed',
              default=False,
              help='Input fastq file is already trimmed for adapter sequence. '
                   'Default: --not-trimmed')
@click.option('--methtype', type=click.STRING,
              default='CG',
              help='Type of methylation that the Percent Methylation bed files '
                   'contain. Choices are: C, CG, CH, CHG, or CHH. Default: CG')
@click.option('--strand', type=click.STRING,
              default='both',
              help='Strand that the Percent Methylation bed files contain '
                   'information about. Choices are: positive, negative, or '
                   'both. Default: both')
@click.option('--max_dup_reads', type=click.INT,
              default=1,
              help='Maximum number of duplicate reads allowed to inform each '
                   'C in the Percent Methylation bed files. Default: 1')
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
              default='info.yaml',
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. '
                   'Default: info.yaml')
@click.argument('input', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def align(in_fastq, out_prefix, out_dir, genome, noadap_bs2_params,
          adaptrim_bs2_params, trimmed, methtype, strand, max_dup_reads,
          threads, working_dir, infoyaml):
    """
    Aligns FASTQ to create BAM permeth BED files.

    Main pipeline for converting raw FASTQ files (from Illumina sequenced
    WGBS samples) into SAM and Percentage Methylation BED format (PerMeth).

    \b
    Required arguments:
    INPUT        Input FASTQ file which will be processed by the pipeline
    OUT_PREFIX   Prefix of all output files
                 This should be short string, not a full path.
                      Ex:  test
                      NOT: test/test
                 If you want to output in a directory other than the current
                 working directory, use Option: out_dir.
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
    #Create output directory
    if out_dir != '':
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    bed_dir = os.path.join(out_dir, 'permethbed_{0}'.format(out_prefix))
    if not os.path.exists(bed_dir):
        os.makedirs(bed_dir)
    bed_prefix = os.path.join(bed_dir, '{}_'.format(out_prefix))

    #Filter fastq file
    logging.info('Filtering out quality failed reads from fastq file')
    fastqtools.qual_filter_fastq(in_fastq, qualfil_fastq)

    #Remove adapter contamination from reads
    logging.info('Removing adapter contamination and split fastq files')
    fastqtools.adapter_remove(qualfil_fastq, noadap_fq, adaptrim_fq,
                             info_dict['adapter'])

    #Align to genome
    logging.info('Aligning reads without adapter contamination to {}'
                 .format(genome))
    noadap_bs2_params = '--bt-p {} {}'.format(threads, noadap_bs2_params)
    bsseeker.align_bs2(bs2_path, noadap_bs2_params, fasta, bs2_index, noadap_fq,
                       noadap_bam)
    logging.info('Aligning reads that had adapter contamination trimmed out '
                 'to {}'.format(genome))
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

    #Index full bam file
    command = 'samtools index {}'.format(full_bam)
    logging.info(command)
    subprocess.check_call(command, shell=True)

    #Convert sam to permeth bed files (percent methylation bed files)
    #This also only prints CpGs on chromosomes in ranges defined in the yaml
    samutils.bam_to_permeth(full_bam, out_prefix, bed_prefix, genome,
                            methtype, strand, max_dup_reads, chroms, threads)

    #Determine conversion efficiency of the experiment


@cli.command()
@click.option('--mask', type=click.STRING,
              default='',
              help='GTF or BED file indicating regions to be masked out of analysis. '
                   'Default: Does not mask any regions.')
@click.option('--raw-data', 'raw_data', type=click.STRING,
              default='',
              help='File name to print raw data for each ROI. This includes '
                   'counts of methylated reads, total reads, and CpGs covered '
                   'in ROI. Default: <Does not produce this file>')
@click.option('--threads', type=click.INT,
              default=NUM_CPUS,
              help='Number of threads used when multiprocessing. '
                   'Default: Number of system CPUs')
@click.option('--min-reads', 'min_read_count', type=click.INT,
              default=1,
              help="Minimum read count for a sample over a given region of "
                   "interest. If this threshold is not met, NA is reported "
                   "for the given sample's ROI. Default: 1")
@click.option('--min-cpgs', 'min_cpg_count', type=click.INT,
              default=10,
              help="Minimum CpG count for a sample over a given region of "
                   "interest. If this threshold is not met, NA is reported "
                   "for the given sample's ROI. Default: 10")
@click.option('--min-samples', 'min_sample_coverage', type=click.INT,
              default=-1,
              help="Minimum number of samples required to report a ROI. For "
                   "example, if this is set to the number of samples input "
                   "and at least one of those samples does not meet the "
                   "minimum read count for the ROI, that ROI is not reported. "
                   "Default: <all samples>")
@click.argument('input_tsv', type=click.STRING)
@click.argument('out_table', type=click.STRING)
@click.argument('roi_file', type=click.STRING)
def roi(input_tsv, out_table, roi_file, mask, min_read_count, min_cpg_count,
        min_sample_coverage, raw_data, threads):
    """
    Calls methylation over ROIs.

    Creates a table that contains the methylation percentage over all Regions of
    Interest (ROIs). Allows creation of an optional 2col table which contains
    the number of methylated and total reads per sample per ROI.

    \b
    Required arguments:
    INPUT_TSV    Input tab separated file indicating sample names and locations.
                 The first tab should be sample name
                 The second tab should be the path to the prefix of bed file
                 Ex file format:
                 pm01   tests/data/bed/pm01_
                 pm02   tests/data/bed/pm02_
    OUT_TABLE    Name of the table which contains all of the ROI methylation
                 information
    ROI_FILE     GTF or BED file indicating the ROI (Regions of Interest).
                 Each ROI will be output as a single line in the OUT_TABLE.
    """
    in_bed_prefixes = []
    in_sample_list = []
    with open(input_tsv, 'r') as infile:
        for line in infile:
            line = line[:-1]
            in_bed_prefixes.append(line.split('\t')[1])
            in_sample_list.append(line.split('\t')[0])
    if min_sample_coverage == -1:
        min_sample_coverage = len(in_sample_list)
    permethbed.roi_meth(in_bed_prefixes, in_sample_list, out_table,
                        mask, roi_file, min_read_count, min_cpg_count,
                        min_sample_coverage, raw_data, threads)


@cli.command()
@click.option('--suffix', type=click.STRING,
              default='',
              help='Requires all in files end in a particular string. Useful '
                   'if you have multiple file types in a folder but only want '
                   'to adjust one kind at a time. Default: None')
@click.option('--col', type=click.INT,
              default=2,
              help='Column number to be changed (0-based). Ex. If set to 2, '
                   'the third column of a file will be changed. If this is a '
                   'bed file, the end location will be changed. Default: 2')
@click.option('--adjust', type=click.INT,
              default=1,
              help='Amount to adjust each number in the column. Ex: If set to '
                   '1, each number in the column will increase by 1. If set '
                   'to -5000, each number will be subtracted by 5000. '
                   'Default: 1')
@click.option('--header/--no-header',
              default=True,
              help='Boolean which indicates if there is a header in the input '
                   'files. Default: --header')
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def adjustcol(in_prefix, out_prefix, suffix, col, adjust, header):
    """
    Adjusts numerical column of files.

    Takes all files matching the prefix of in_prefix and outputs

    \b
    Required arguments:
    IN_PREFIX    Prefix of all files input into the
    OUT_PREFIX   Prefix of all output files
    """
    for in_file_name in glob.glob('{}*{}'.format(in_prefix, suffix)):
        suffix = in_file_name.split(in_prefix)[1]
        out_file_name = '{}{}'.format(out_prefix, suffix)
        print('Processing: {}'.format(in_file_name))
        print('Output to:  {}'.format(out_file_name))
        outfile = open(out_file_name, 'wb')
        headerin = header
        with open(in_file_name, 'r') as in_file:
            for line in in_file:
                if headerin:
                    outfile.write(line)
                    headerin = False
                else:
                    line = line[:-1]
                    cells = line.split('\t')
                    if len(cells) > col:
                        cells[col] = int(cells[col]) + adjust
                        printline = str(cells.pop(0))
                        for cell in cells:
                            printline = '{}\t{}'.format(printline, cell)
                        printline = '{}\n'.format(printline)
                        outfile.write(printline)
                    else:
                        print('Warning, col {} does not exist in line:\n{}'
                              .format(col, line))
        outfile.close()


@cli.command()
@click.option('--window-size', 'windowsize', type=click.INT,
              default=20000,
              help='Size of windows in bp. Default: 20000')
@click.option('--mask', type=click.STRING,
              default='',
              help='GTF or BED file indicating regions to be masked out of analysis. '
                   'Default: Does not mask any regions.')
@click.option('--raw-data', 'raw_data', type=click.STRING,
              default='',
              help='File name to print raw data for each window. This includes '
                   'counts of methylated reads, total reads, and CpGs covered '
                   'in each window. Default: <Does not produce this file>')
@click.option('--threads', type=click.INT,
              default=NUM_CPUS,
              help='Number of threads used when multiprocessing. '
                   'Default: Number of system CPUs')
@click.option('--min-reads', 'min_read_count', type=click.INT,
              default=1,
              help="Minimum read count for a sample over a given region of "
                   "interest. If this threshold is not met, NA is reported "
                   "for the given sample's window. Default: 1")
@click.option('--min-cpgs', 'min_cpg_count', type=click.INT,
              default=20,
              help="Minimum CpG count for a sample over a given region of "
                   "interest. If this threshold is not met, NA is reported "
                   "for the given sample's window. Default: 20")
@click.option('--min-samples', 'min_sample_coverage', type=click.INT,
              default=-1,
              help="Minimum number of samples required to report a window. For "
                   "example, if this is set to the number of samples input "
                   "and at least one of those samples does not meet the "
                   "minimum read count for the window, that window is not "
                   "reported. Default: <all samples>")
@click.option('--genome', type=click.STRING,
              default='hg38',
              help='Genome used for alignment and analysis. '
                   'Default: hg38')
@click.option('--infoyaml', type=click.STRING,
              default='info.yaml',
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. '
                   'Default: info.yaml')
@click.argument('input_tsv', type=click.STRING)
@click.argument('out_table', type=click.STRING)
def window(input_tsv, out_table, windowsize, mask, raw_data, threads,
           min_read_count, min_cpg_count, min_sample_coverage, genome,
           infoyaml):
    """
    Calls methylation over windows.

    Creates a table that contains windows of methylation percentage over
    all areas/chromosomes of the genome. Allows creation of an optional 2col
    table which contains the number of methylated and total reads per sample
    per window.

    \b
    Required arguments:
    INPUT_TSV    Input tab separated file indicating sample names and locations.
                 The first tab should be sample name
                 The second tab should be the path to the prefix of bed file
                 Ex file format:
                 pm01   tests/data/bed/pm01_
                 pm02   tests/data/bed/pm02_
    OUT_TABLE    Name of the table which contains all of the window methylation
                 information
    """
    workingdir = tempfile.mkdtemp()
    window_roi = os.path.join(workingdir, 'windows.bed')
    in_bed_prefixes = []
    in_sample_list = []
    with open(input_tsv, 'r') as infile:
        for line in infile:
            line = line[:-1]
            in_bed_prefixes.append(line.split('\t')[1])
            in_sample_list.append(line.split('\t')[0])
    if min_sample_coverage == -1:
        min_sample_coverage = len(in_sample_list)
    stream = file(infoyaml, 'r')
    info_dict = yaml.safe_load(stream)
    chroms = info_dict[genome]['chroms']
    permethbed.create_window_roi(window_roi, windowsize, chroms)
    permethbed.roi_meth(in_bed_prefixes, in_sample_list, out_table,
                        mask, window_roi, min_read_count, min_cpg_count,
                        min_sample_coverage, raw_data, threads)


@cli.command()
@click.argument('input_pm_prefix', type=click.STRING)
@click.argument('out_dss_prefix', type=click.STRING)
def pm2dss(in_prefix, out_prefix):
    """"""
    bed_files = []
    for file_name in glob.glob('{}*.bed'.format(in_prefix)):
        bed_files.append(file_name)
    gz_bed_files = []
    for file_name in glob.glob('{}*.bed.gz'.format(in_prefix)):
        gz_bed_files.append(file_name)

