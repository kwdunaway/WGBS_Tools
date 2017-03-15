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
from pkg_resources import resource_filename
import wgbs_tools
from pybedtools import BedTool

from wgbs_tools import fastqtools
from wgbs_tools import bsseeker
from wgbs_tools import utilities
from wgbs_tools import samutils
from wgbs_tools import permethbed

NUM_CPUS = multiprocessing.cpu_count()

@click.group()
def cli():
    """A versatile toolkit to manipulate and analyze WGBS data."""
    pass


@cli.command()
@click.option('--out-dir', 'out_dir', type=click.STRING,
              default='',
              help='Directory to put all outfiles. '
                   'Default: <directory with same name as OUT_PREFIX>')
@click.option('--genome', type=click.STRING,
              default='hg38',
              help='Genome used for alignment and analysis. '
                   'Default: hg38')
@click.option('--chew', 'chew_length', type=click.INT,
              default=10,
              help='Number of bases to removed off of the end of each read. '
                   'Default: 10')
@click.option('--min-read', 'min_seqlength', type=click.INT,
              default=35,
              help='Minimum read length of reads after adapter trimming and '
                   'chew. If the read is less than this lenght, it is not '
                   'included in the output. Default: 35')
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
@click.option('--conv-chr', 'conv_chrom', type=click.STRING,
              default='chrM',
              help='Chromosome that determines conversion efficency. Default '
                   'is the mitochondria chromosome, although if unmethylated '
                   'DNA (ex: lamda) was spiked in your sample before bisulfite '
                   'conversion, then you should use that. Default: chrM')
@click.option('--threads', type=click.INT,
              default=NUM_CPUS,
              help='Number of threads used when multiprocessing. '
                   'Default: Number of system CPUs')
@click.option('--working-dir', 'working_dir', type=click.STRING,
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
@click.argument('in_fastq', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def process_se(in_fastq, out_prefix, out_dir, chew_length, min_seqlength,
               genome, noadap_bs2_params, adaptrim_bs2_params, methtype,
               strand, max_dup_reads, conv_chrom, threads, working_dir,
               infoyaml):
    """
    Pipeline to process single end FASTQ file.

    Main pipeline for converting a single end sequenced raw FASTQ files (from
    Illumina sequenced WGBS samples) into BAM and Percentage Methylation BED
    format (PerMeth).

    \b
    Required arguments:
    IN_FASTQ        Input FASTQ file which will be processed by the pipeline
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
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        workingdir = working_dir
    #Load data
    if infoyaml == 'info.yaml':
        infoyaml = resource_filename(wgbs_tools.__name__, '../info.yaml')
    stream = file(infoyaml, 'r')
    info_dict = yaml.safe_load(stream)
    bs2_path = info_dict['bs2_path']
    index = info_dict[genome]['index']
    fasta = info_dict[genome]['fasta']
    chroms = info_dict[genome]['chroms']
    adapter = info_dict['adapter']

    #Create output directory and file names that go there
    if out_dir == '':
        out_dir = out_prefix
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    full_bam = os.path.join(out_dir, '{}.bam'.format(out_prefix))
    bed_dir = os.path.join(out_dir, 'permethbed_{}'.format(out_prefix))
    if not os.path.exists(bed_dir):
        os.makedirs(bed_dir)
    bed_prefix = os.path.join(bed_dir, '{}_'.format(out_prefix))
    conv_eff = os.path.join(out_dir, '{}_conveff.txt'.format(out_prefix))
    out_summary = os.path.join(out_dir, '{}_summary.txt'.format(out_prefix))

    #Name temp files
    temp_prefix = os.path.join(workingdir, out_prefix)
    qualfil_fastq = '{}_filtered.fq.gz'.format(temp_prefix)
    noadap_fq = '{}_noadap.fq.gz'.format(temp_prefix)
    adaptrim_fq = '{}_trimmed.fq.gz'.format(temp_prefix)
    noadap_bam = '{}_noadap.bam'.format(temp_prefix)
    noadap_log = '{}_noadap.bam.bs_seeker2_log'.format(temp_prefix)
    adaptrim_bam = '{}_adaptrim.bam'.format(temp_prefix)
    adaptrim_log = '{}_adaptrim.bam.bs_seeker2_log'.format(temp_prefix)
    noadap_sorted = '{}_noadap_sorted'.format(temp_prefix)
    adaptrim_sorted = '{}_adaptrim_sorted'.format(temp_prefix)

    #Filter fastq file
    logging.info('Filtering out quality failed reads from fastq file')
    fastqtools.qual_filter_fastq(in_fastq, qualfil_fastq)

    #Remove adapter contamination from reads
    logging.info('Removing adapter contamination and split fastq files')
    fastqtools.adapter_remove(qualfil_fastq, noadap_fq, adaptrim_fq, adapter,
                              chew_length, min_seqlength)

    #Align to genome
    logging.info('Aligning reads without adapter contamination to {}'
                 .format(genome))
    noadap_bs2_params = '--bt-p {} {}'.format(threads, noadap_bs2_params)
    bsseeker.align_bs2(bs2_path, noadap_bs2_params, fasta, index, noadap_fq,
                       noadap_bam)
    logging.info('Aligning reads that had adapter contamination trimmed out '
                 'to {}'.format(genome))
    adaptrim_bs2_params = '--bt-p {} {}'.format(threads, adaptrim_bs2_params)
    bsseeker.align_bs2(bs2_path, adaptrim_bs2_params, fasta, index,
                       adaptrim_fq, adaptrim_bam)

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

    #Determine conversion efficiency of the experiment if designated
    #     chromosome is in chroms list.
    if conv_chrom in chroms:
        conv_bed = '{}_chrM.bed.gz'.format(bed_prefix)
        conv_eff_dict = permethbed.bed_meth_stats(conv_bed)
        ce = open(conv_eff, 'wb')
        header_line = 'conv_eff\tmeth\ttotal\tcpg_count\n'
        ce.write(header_line)
        eff = 1 - conv_eff_dict['perc']
        info_line = '{}\t{}\t{}\t{}\n'\
            .format(eff, conv_eff_dict['meth'], conv_eff_dict['total'],
                    conv_eff_dict['cpg'])
        ce.write(info_line)
    else:
        logging.warning('Conversion efficency was not calculated because '
                        '{} is not one of the designated '
                        'chromosomes.'.format(conv_chrom))
        conv_eff = ''

    #Create summary file from alignment logs and conversion efficency
    bsseeker.process_logs(noadap_log, adaptrim_log, conv_eff, out_summary)

    #Remove temp files
    if working_dir == '':
        shutil.rmtree(workingdir)


@cli.command()
@click.option('--out_dir', type=click.STRING,
              default='',
              help='Directory to put all outfiles. '
                   'Default: <current working directory>')
@click.option('--genome', type=click.STRING,
              default='hg38',
              help='Genome used for alignment and analysis. '
                   'Default: hg38')
@click.option('--chew', 'chew_length', type=click.INT,
              default=10,
              help="Length in bp that the 3' end of a read is removed. "
                   "Default: 10")
@click.option('--read-min', 'min_readlength', type=click.INT,
              default=35,
              help="Minimum read length allowed after adapter trimming and "
                   "chew. Any reads shorter will be discarded. Default: 35")
@click.option('--max_dup_reads', type=click.INT,
              default=1,
              help='Maximum number of duplicate reads allowed to inform each '
                   'C in the Percent Methylation bed files. Default: 1')
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
@click.option('--conv-chr', 'conv_chrom', type=click.STRING,
              default='chrM',
              help='Chromosome that determines conversion efficency. Default '
                   'is the mitochondria chromosome, although if unmethylated '
                   'DNA (ex: lamda) was spiked in your sample before bisulfite '
                   'conversion, then you should use that. Default: chrM')
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
@click.argument('in_fastq_f', type=click.STRING)
@click.argument('in_fastq_r', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def process_pe(in_fastq_f, in_fastq_r, out_prefix, out_dir, genome, chew_length,
               min_readlength, noadap_bs2_params, adaptrim_bs2_params,
               methtype, strand, max_dup_reads, conv_chrom, threads,
               working_dir, infoyaml):
    """
    Pipeline to process paired end FASTQ files.

    Main pipeline for converting a paired end FASTQ experiment (from Illumina
    sequenced WGBS samples) into BAM and Percentage Methylation BED format
    (PerMeth).

    \b
    Required arguments:
    IN_FASTQ_F   Input forward FASTQ file
    IN_FASTQ_R   Input reverse FASTQ file
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
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        workingdir = working_dir
    #Load data
    if infoyaml == 'info.yaml':
        infoyaml = resource_filename(wgbs_tools.__name__, '../info.yaml')
    stream = file(infoyaml, 'r')
    info_dict = yaml.safe_load(stream)
    bs2_path = info_dict['bs2_path']
    index = info_dict[genome]['index']
    fasta = info_dict[genome]['fasta']
    chroms = info_dict[genome]['chroms']
    adapter = info_dict['adapter']

    #Name temp files
    temp_prefix = os.path.join(workingdir, out_prefix)
    noadap_fq_f = '{}_noadap_f.fq.gz'.format(temp_prefix)
    adaptrim_fq_f = '{}_trimmed_f.fq.gz'.format(temp_prefix)
    noadap_fq_r = '{}_noadap_r.fq.gz'.format(temp_prefix)
    adaptrim_fq_r = '{}_trimmed_r.fq.gz'.format(temp_prefix)
    noadap_bam = '{}_noadap.bam'.format(temp_prefix)
    noadap_log = '{}_noadap.bam.log'.format(temp_prefix)
    adaptrim_bam = '{}_adaptrim.bam'.format(temp_prefix)
    adaptrim_log = '{}_adaptrim.bam.log'.format(temp_prefix)
    noadap_sorted = '{}_noadap_sorted'.format(temp_prefix)
    adaptrim_sorted = '{}_adaptrim_sorted'.format(temp_prefix)

    #Create output directory and file names that go there
    if out_dir != '':
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    full_bam = os.path.join(out_dir, '{}.bam'.format(out_prefix))
    bed_dir = os.path.join(out_dir, 'permethbed_{}'.format(out_prefix))
    if not os.path.exists(bed_dir):
        os.makedirs(bed_dir)
    bed_prefix = os.path.join(bed_dir, '{}_'.format(out_prefix))
    conv_eff = os.path.join(out_dir, '{}_conveff.txt'.format(out_prefix))
    out_summary = os.path.join(out_dir, '{}_summary.txt'.format(out_prefix))

    #Remove adapter contamination from reads
    logging.info('Removing adapter contamination and split fastq files')
    fastqtools.pe_adapter_remove(in_fastq_f, noadap_fq_f, adaptrim_fq_f,
                                 adapter, in_fastq_r, noadap_fq_r, adaptrim_fq_r,
                                 adapter, chew_length, min_readlength, threads)

    #Align to genome
    logging.info('Aligning reads without adapter contamination to {}'
                 .format(genome))
    noadap_bs2_params = '--bt-p {} {}'.format(threads, noadap_bs2_params)
    bsseeker.align_bs2(bs2_path, noadap_bs2_params, fasta, index, noadap_fq,
                       noadap_bam)
    logging.info('Aligning reads that had adapter contamination trimmed out '
                 'to {}'.format(genome))
    adaptrim_bs2_params = '--bt-p {} {}'.format(threads, adaptrim_bs2_params)
    bsseeker.align_bs2(bs2_path, adaptrim_bs2_params, fasta, index,
                       adaptrim_fq, adaptrim_bam)

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

    #Determine conversion efficiency of the experiment if designated
    #     chromosome is in chroms list.
    if conv_chrom in chroms:
        conv_bed = '{}_chrM.bed.gz'.format(bed_prefix)
        conv_eff_dict = permethbed.bed_meth_stats(conv_bed)
        ce = open(conv_eff, 'wb')
        header_line = 'conv_eff\tmeth\ttotal\tcpg_count\n'
        ce.write(header_line)
        eff = 1 - conv_eff_dict['perc']
        info_line = '{}\t{}\t{}\t{}\n'\
            .format(eff, conv_eff_dict['meth'], conv_eff_dict['total'],
                    conv_eff_dict['cpg'])
        ce.write(info_line)
    else:
        logging.warning('Conversion efficency was not calculated because '
                        '{} is not one of the designated '
                        'chromosomes.'.format(conv_chrom))
        conv_eff = ''

    #Create summary file from alignment logs and conversion efficency
    bsseeker.process_logs(noadap_log, adaptrim_log, conv_eff, out_summary)

    #Remove temp files
    if working_dir != '':
        shutil.rmtree(workingdir)


@cli.command()
@click.option('--conv', 'conv_eff', type=click.STRING,
              default='',
              help='Conversion efficency summary file. If none entered then '
                   'those columns will not be included in the output summary '
                   'file. Default: <None>')
@click.argument('noadap_log', type=click.STRING)
@click.argument('adaptrim_log', type=click.STRING)
@click.argument('out_summary', type=click.STRING)
def sumlogs(noadap_log, adaptrim_log, out_summary, conv_eff):
    """
    Summarizes BS Seeker2 logs.

    Takes in both a noadap_log file and an adaptrim_log file. Then, parses
    the relevant information and outputs it into a simple tab separated
    table. This also allows conversion efficency to be added to the rightmost
    columns.

    \b
    Required arguments:
    NOADAP_LOG     Log file for the noadapter BS Seeker2 alignment
    ADAPTRIM_LOG   Log file for the adapter trimmed BS Seeker2 alignment
    OUT_SUMMARY    Name of the output file
     """
    bsseeker.process_logs(noadap_log, adaptrim_log, conv_eff, out_summary)


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
@click.option('--debug/--no-debug',
              default=False,
              help='Print debug messages. Default: --no-debug')
@click.argument('input_tsv', type=click.STRING)
@click.argument('out_table', type=click.STRING)
@click.argument('roi_file', type=click.STRING)
def roi(input_tsv, out_table, roi_file, mask, min_read_count, min_cpg_count,
        min_sample_coverage, raw_data, threads, debug):
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
    if debug:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.DEBUG)
    in_bed_prefixes = []
    in_sample_list = []
    with open(input_tsv, 'r') as infile:
        for line in infile:
            # if line.endswith('\n') or line.endswith('\r'):
            line = line.rstrip()
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
@click.option('--cols', type=click.STRING,
              default='2',
              help='Column number(s) to be changed (0-based). Ex. If set to 2, '
                   'the third column of a file will be changed. If this is set '
                   'to 1,2 then the second and third columns are changed. '
                   'Default: 2')
@click.option('--adjusts', type=click.STRING,
              default='1',
              help='Amount to adjust each number in the column. Ex: If set to '
                   '1, each number in the column will increase by 1. If set '
                   'to 1,-5000, each number in the first designated column '
                   'will increase by 1 and the second designated column will '
                   'have 5000 subtracted from it. Default: 1')
@click.option('--header/--no-header',
              default=True,
              help='Boolean which indicates if there is a header in the input '
                   'files. Default: --header')
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def adjustcols(in_prefix, out_prefix, suffix, cols, adjusts, header):
    """
    Adjusts numerical column of files.

    Takes all files matching the prefix of in_prefix and outputs

    \b
    Required arguments:
    IN_PREFIX    Prefix of all files input into the
    OUT_PREFIX   Prefix of all output files
    """
    col_list = cols.split(',')
    adjust_list = adjusts.split(',')
    assert len(col_list) == len(adjust_list), \
        "Error: You are trying to adjust {} columns with {} adjustments"\
            .format(len(col_list), len(adjust_list))
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
                    line = line.rstrip()
                    cells = line.split('\t')
                    for num in range(len(col_list)):
                        col = int(col_list[num])
                        adjust = int(adjust_list[num])
                        if len(cells) > col:
                            cells[col] = int(cells[col]) + adjust
                        else:
                            print('Warning, col {} does not exist in line:\n{}'
                                  .format(col, line))
                    printline = str(cells.pop(0))
                    for cell in cells:
                        printline = '{}\t{}'.format(printline, cell)
                    printline = '{}\n'.format(printline)
                    outfile.write(printline)
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
            line = line.rstrip()
            line_list = line.split('\t')
            if len(line_list) > 1:
                in_bed_prefixes.append(line_list[1])
                in_sample_list.append(line_list[0])
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
@click.option('--gz/--no-gz',
              default=True,
              help='Boolean which indicates if the output files will be '
                   'compressed (.gz format) or not. Default: --gz (compressed)')
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def pm2dss(in_prefix, out_prefix, gz):
    """
    Converts pm_bed to dss format.

    The percent methylation bed files (pm_bed) are used to create DSS files,
    which are used in multiple R packages, including DSSfinder. The output
    can be either regular or gzipped files (for space conservation).

    \b
    Required arguments:
    IN_PREFIX      Prefix for all percent methlyated bed files (can be
                   compressed and/or uncompressed)
    OUT_PREFIX     Prefix for all outfiles (in DSS format). The full name of
                   each outfile is [out_prefix][uniq_id].dss

    \b
    Example run: wgbs_tools pm2dss pm01_bed/pm01 pm01_dss/pm01_
    If there were three files in the folder:
        pm01_bed/pm01_chr1.bed
        pm01_bed/pm01_chr2.bed
        pm01_bed/pm01_chr3.bed.gz
    The output would be put in these three files:
        pm01_dss/pm01_chr1.dss
        pm01_dss/pm01_chr2.dss
        pm01_dss/pm01_chr3.dss
    """
    if gz:
        suffix = '.dss.gz'
    else:
        suffix = '.dss'
    for bed_name in glob.glob('{}*.bed'.format(in_prefix)):
        uniqname = bed_name.split(in_prefix)[1].split('.bed')[0]
        dss_name = '{}{}{}'.format(out_prefix, uniqname, suffix)
        logging.info('Processing {}'.format(bed_name))
        logging.info('Creating {}'.format(dss_name))
        permethbed.convert_pm2dss(bed_name, dss_name)
    for bed_name in glob.glob('{}*.bed.gz'.format(in_prefix)):
        uniqname = bed_name.split(in_prefix)[1].split('.bed.gz')[0]
        dss_name = '{}{}{}'.format(out_prefix, uniqname, suffix)
        logging.info('Processing {}'.format(bed_name))
        logging.info('Creating {}'.format(dss_name))
        permethbed.convert_pm2dss(bed_name, dss_name)


@cli.command()
# @click.option('--gz/--no-gz',
#               default=True,
#               help='Boolean which indicates if the output files will be '
#                    'compressed (.gz format) or not. Default: --gz (compressed)')
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def pm2bg(in_prefix, out_file, gz):
    """
    Converts pm_bed to bedgraph format.

    The percent methylation bed files (pm_bed) are used to create a bedGraph
    file, which can be used for various reasons, including creating a trackhub.

    \b
    Required arguments:
    IN_PREFIX      Prefix for all percent methlyated bed files (can be
                   compressed and/or uncompressed)
    OUT_FILE       Prefix for all outfiles (in bedgraph format). The full
                   name of each outfile is [out_prefix][uniq_id].bg

    \b
    Example run: wgbs_tools pm2bg pm01_bed/pm01 pm01_bg/pm01_
    If there were three files in the folder:
        pm01_bed/pm01_chr1.bed
        pm01_bed/pm01_chr3.bed.gz
    The output would be put in these three files:
        pm01_bg/pm01_chr1.bg
        pm01_bg/pm01_chr3.bg
    """
    #TODO: convert this to take in multiple permeth files and output a single bg
    # if gz:
    #     suffix = '.bg.gz'
    # else:
    #     suffix = '.bg'
    # for bed_name in glob.glob('{}*.bed'.format(in_prefix)):
    #     uniqname = bed_name.split(in_prefix)[1].split('.bed')[0]
    #     bg_name = '{}{}{}'.format(out_prefix, uniqname, suffix)
    #     logging.info('Processing {}'.format(bed_name))
    #     logging.info('Creating {}'.format(bg_name))
    #     permethbed.convert_pm2bg(bed_name, bg_name)
    # for bed_name in glob.glob('{}*.bed.gz'.format(in_prefix)):
    #     uniqname = bed_name.split(in_prefix)[1].split('.bed.gz')[0]
    #     bg_name = '{}{}{}'.format(out_prefix, uniqname, suffix)
    #     logging.info('Processing {}'.format(bed_name))
    #     logging.info('Creating {}'.format(bg_name))
    #     permethbed.convert_pm2bg(bed_name, bg_name)
    #

@cli.command()
@click.option('--suffix', type=click.STRING,
              default='',
              help='Requires all in files end in a particular string. Useful '
                   'if you have multiple file types in a folder but only want '
                   'to adjust one kind at a time. Default: None')
@click.argument('in_prefix', type=click.STRING)
def ll_fixdmrs(in_prefix, suffix):
    """
    Fixes all single line DMR files in a folder.

    Takes all files matching the prefix of in_prefix and searches for lines
    with a single column in the second line. If that is the case, it takes
    all of the lines in the file (after the header) and combines them into a
    single line that is tab separated. This overwrites the previous file.

    \b
    Required arguments:
    IN_PREFIX    Prefix of all files input into the
    """
    workingdir = tempfile.mkdtemp()
    for in_file_name in glob.glob('{}*{}'.format(in_prefix, suffix)):
        suffix = in_file_name.split(in_prefix)[1]
        out_file_name = os.path.join(workingdir, 'temp_{}'.format(suffix))
        processed = 'no'
        with open(in_file_name, 'r') as in_file:
            try:
                logging.warning('{} does not have header line. Skipping...'
                                .format(in_file_name))
                header = next(in_file)
            except:
                continue
            try:
                logging.warning('{} does not have a second line. Skipping...'
                                .format(in_file_name))
                firstline = next(in_file)
            except:
                continue
            if len(firstline.split('\t')) == 1:
                print('Processing: {}'.format(in_file_name))
                firstline = firstline[:-1]
                for line in in_file:
                    line = line.rstrip()
                    firstline = '{}\t{}'.format(firstline, line)
                firstline = '{}\n'.format(firstline)
                outfile = open(out_file_name, 'wb')
                outfile.write(header)
                outfile.write(firstline)
                outfile.close()
                processed = 'yes'
        if processed == 'yes':
            command = 'rm {}'.format(in_file_name)
            subprocess.check_call(command, shell=True)
            command = 'mv {} {}'.format(out_file_name, in_file_name)
            subprocess.check_call(command, shell=True)


@cli.command()
@click.option('--suffix', type=click.STRING,
              default='.bed.gz',
              help='Requires all in files end in a particular string. '
                   'Default: .bed.gz')
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_file', type=click.STRING)
def pm_stats(in_prefix, out_file, suffix):
    """
    Gets stats of multiple pm_bed files.

    \b
    Prints 5 columns:
    [0] name: Unique name of file (string between prefix and suffix)
    [1] percentage: Percentage methylation of file
    [2] methylated_reads: Amount of methylated reads in pm file
    [3] total_reads: Amount of total reads in pm file
    [4] cpg_count: Number of CpGs the pm file has information for

    \b
    This can also be used to determine conversion efficiency by running it on
    the conv_eff chromosome (see below for explanation). Then subtract the
    percentage from 1 (ie: 1 - percentage). There are 3 possible conv_eff
    chromosomes:
      1.  If you spiked your samples with lamda DNA (known to be 100%
          unmethylated) prior to bisulfite conversion, run this on that
          chromosome. Note, you will need to ensure you added the lamda
          sequence to info.yaml and the BS_Seeker index prior to alignment.
      2.  chrM. This is the most likely choice if you didn't spike your samples.
      3.  CH pm_bed files. Create pm_bed files looking at CH methylation and
          then run this command on all of those files. This assumes there is
          no CH methylation (which is known to exist in plants and certain
          mammalian tissues).

    \b
    Required arguments:
    IN_PREFIX      Prefix for all percent methlyated bed files (can be
                   compressed and/or uncompressed)
    OUT_FILE       Tab separated table containing statistical information
    """
    out = open(out_file, 'wb')
    outline = 'name\tpercentage\tmethylated_reads\ttotal_reads\tcpg_count\n'
    out.write(outline)
    for in_file_name in glob.glob('{}*{}'.format(in_prefix, suffix)):
        uniqname = in_file_name.split(in_prefix)[1].split(suffix)[0]
        meth_dict = permethbed.bed_meth_stats(in_file_name)
        outline = '{}\t{}\t{}\t{}\t{}\n' \
            .format(uniqname, meth_dict['perc'], meth_dict['meth'],
                    meth_dict['total'], meth_dict['cpg'])
        out.write(outline)
    out.close()


@cli.command()
@click.option('--adapter', type=click.STRING,
              default='AGATCGGAAG',
              help='Beginning sequence of adapter. Default: AGATCGGAAG')
@click.option('--out_dir', type=click.STRING,
              default='',
              help='Directory to put all outfiles. '
                   'Default: <current working directory>')
@click.option('--threads', type=click.INT,
              default=NUM_CPUS,
              help='Number of threads used when multiprocessing. '
                   'Default: Number of system CPUs')
@click.option('--chew', 'chew_length', type=click.INT,
              default=10,
              help='Number of bases to removed off of the end of each read. '
                   'Default: 10')
@click.option('--min-read', 'min_seqlength', type=click.INT,
              default=35,
              help='Minimum read length of reads after adapter trimming and '
                   'chew. If the read is less than this lenght, it is not '
                   'included in the output. Default: 35')
@click.option('--working_dir', type=click.STRING,
              default='',
              help='Working directory where temp files are written and read. '
                   'Default: <EMPTY> (Uses tempfile.mkdtemp() to define a '
                   'temperary working directory)')
@click.argument('in_fastq', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def trim_sefq(in_fastq, out_prefix, adapter, out_dir, threads, chew_length,
              min_seqlength, working_dir):
    """
    Filters and trims single end FASTQ file.

    \b
    Quality filters and adapter trims a pair of paired-end FASTQ files. Takes in
    two FASTQ file and outputs four different FASTQ files:
      1) *_trimmed.fq.gz Contains reads that had adapter sequence detected and
                         have been trimmed out. Then, the sequence was chewed
                         back another 10 bp.
      2) *_noadap.fq.gz  Contains reads that had no adapter sequence detected.
    \b
    Required arguments:
    IN_FASTQ         Input FASTQ file
    OUT_PREFIX       Prefix of the two output files
    """
    #Create output directory
    if out_dir != '':
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    if working_dir == '':
        workingdir = tempfile.mkdtemp()
    else:
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        workingdir = working_dir

    #Name files
    temp_prefix = out_prefix.split('\t')[-1]
    qualfil_fastq = '{}_filtered.fq.gz'.format(temp_prefix)
    noadap_fq = '{}_noadap.fq.gz'.format(temp_prefix)
    adaptrim_fq = '{}_trimmed.fq.gz'.format(temp_prefix)

    #Filter fastq file
    logging.info('Filtering out quality failed reads from fastq file')
    fastqtools.qual_filter_fastq(in_fastq, qualfil_fastq)

    #Remove adapter contamination from reads
    logging.info('Removing adapter contamination and split fastq files')
    fastqtools.adapter_remove(qualfil_fastq, noadap_fq, adaptrim_fq, adapter,
                              chew_length, min_seqlength)

    if working_dir == '':
        shutil.rmtree(workingdir)


@cli.command()
@click.option('--for-adap', 'for_adap', type=click.STRING,
              default='AGATCGGAAG',
              help='Beginning sequence of forward adapters. Default: '
                   'AGATCGGAAG')
@click.option('--rev-adap', 'rev_adap', type=click.STRING,
              default='AGATCGGAAG',
              help='Beginning sequence of forward adapters. Default: '
                   'AGATCGGAAG')
@click.option('--out_dir', type=click.STRING,
              default='',
              help='Directory to put all outfiles. '
                   'Default: <current working directory>')
@click.option('--threads', type=click.INT,
              default=NUM_CPUS,
              help='Number of threads used when multiprocessing. '
                   'Default: Number of system CPUs')
@click.option('--chew', 'chew_length', type=click.INT,
              default=10,
              help='Number of bases to removed off of the end of each read. '
                   'Default: 10')
@click.option('--min-read', 'min_seqlength', type=click.INT,
              default=35,
              help='Minimum read length of reads after adapter trimming and '
                   'chew. If the read is less than this lenght, it is not '
                   'included in the output. Default: 35')
@click.argument('in_for_fq', type=click.STRING)
@click.argument('in_rev_fq', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def trim_pefq(in_for_fq, in_rev_fq, out_prefix, for_adap, rev_adap, out_dir,
              threads, chew_length, min_seqlength):
    """
    Filters and trims paired end FASTQ files.

    \b
    Quality filters and adapter trims a pair of paired-end FASTQ files. Takes in
    two FASTQ file and outputs four different FASTQ files:
      1) *_F_trimmed.fq.gz Contains reads that had adapter sequence detected and
                           have been trimmed out. Then, the sequence was chewed
                           back another 10 bp. Forward reads.
      2) *_F_noadap.fq.gz  Contains reads that had no adapter sequence detected.
                           Forward reads.
      3) *_R_trimmed.fq.gz Contains reads that had adapter sequence detected and
                           have been trimmed out. Then, the sequence was chewed
                           back another 10 bp. Reverse reads.
      4) *_R_noadap.fq.gz  Contains reads that had no adapter sequence detected.
                           Reverse reads.

    \b
    Required arguments:
    IN_FOR_FQ        Input Forward orientation FASTQ file
    IN_FOR_FQ        Input Reverse orientation FASTQ file
    OUT_PREFIX       Prefix of the four output files
    """
    #Name files
    fnoadap_fq = '{}_F_noadap.fq.gz'.format(out_prefix)
    fadaptrim_fq = '{}_F_trimmed.fq.gz'.format(out_prefix)
    rnoadap_fq = '{}_R_noadap.fq.gz'.format(out_prefix)
    radaptrim_fq = '{}_R_trimmed.fq.gz'.format(out_prefix)

    #Create output directory
    if out_dir != '':
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    #Remove adapter contamination from reads
    logging.info('Removing adapter contamination and split fastq files')
    fastqtools.pe_adapter_remove(in_for_fq, fnoadap_fq, fadaptrim_fq,
                                 for_adap, in_rev_fq, rnoadap_fq, radaptrim_fq,
                                 rev_adap, chew_length, min_seqlength, threads)


@cli.command()
@click.option('--genome', type=click.STRING,
              default='hg38',
              help='Genome used for alignment and analysis. '
                   'Default: hg38')
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
@click.option('--header', 'header_name', type=click.STRING,
              default='',
              help='Name in header of bed file. This is useful for loading '
                   'these files into a genome browser. The name of the track '
                   'will be {header}_chr#. '
                   'Default: bed_prefix of lowest directory')
@click.option('--infoyaml', type=click.STRING,
              default='info.yaml',
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. '
                   'Default: info.yaml')
@click.argument('in_bam', type=click.STRING)
@click.argument('bed_prefix', type=click.STRING)
def bam2pm(in_bam, bed_prefix, genome, methtype, strand, max_dup_reads, threads,
           header_name, infoyaml):
    """
    Converts BAM to percent methylation BED.

    Converts a BAM file (produced by BS_Seeker2) into a Percentage Methylation
    BED format (PerMeth).

    \b
    Required arguments:
    IN_BAM       Input sorted BAM file.
                 It should be indexed as well, although this will index the
                 bam file if an index is not found.
    BED_PREFIX   Prefix of output percent methylation bed files.
                 If you want to output in a directory other than the current
                 working directory, use Option: out_dir.
    """
    #Load data
    if infoyaml == 'info.yaml':
        infoyaml = resource_filename(wgbs_tools.__name__, '../info.yaml')
    stream = file(infoyaml, 'r')
    info_dict = yaml.safe_load(stream)
    chroms = info_dict[genome]['chroms']

    if not bed_prefix.endswith('_'):
        bed_prefix = '{}_'.format(bed_prefix)

    #Index full bam file if not found
    indexname = '{}.bai'.format(in_bam)
    if not os.path.isfile(indexname):
        command = 'samtools index {}'.format(in_bam)
        logging.warning('Index not found, running command: {}'.format(command))
        subprocess.check_call(command, shell=True)

    #Determine header_prefix
    if header_name == '':
        bed_dirs = bed_prefix.split('/')
        header_prefix = bed_dirs[-1]
    else:
        header_prefix = header_name

    #Convert sam to permeth bed files (percent methylation bed files)
    #This also only prints CpGs on chromosomes in ranges defined in the yaml
    samutils.bam_to_permeth(in_bam, bed_prefix, header_prefix, genome,
                            methtype, strand, max_dup_reads, chroms, threads)


@cli.command()
@click.argument('in_bed', type=click.STRING)
@click.argument('out_table', type=click.STRING)
def ll_chrcov(in_bed, out_table):
    """
    Gets chromosome coverage of bed file.

    Takes in bed file(s) and outputs a table with the chromosome coverage.
    This just yields the total bases covered. If you want a percentage of the
    chromsome covered, you need to divide these numbers by the chromosome
    length.

    \b
    Required arguments:
    IN_BED     Bed file(s) that you want to get chromosome coverage for. If
               you want to input multiple bed files, you must comma separate
               them. For example:
               samp1.bed,samp2.bed,samp3.bed
    OUT_TABLE  File name of the output table. Every line will be a chromosome
               while each column will represent a sample.
    """
    cov_table = {}
    in_beds = in_bed.split(',')
    outfile = open(out_table, 'wb')
    outheader = 'chromosome'
    for bed_file in in_beds:
        cov_table[bed_file] = {}
        bed = BedTool(bed_file)
        outheader = '{}\t{}'.format(outheader, bed_file)
        for feature in bed:
            if feature.chrom in cov_table[bed_file]:
                cov_table[bed_file][feature.chrom] += feature.end - feature.start
            else:
                cov_table[bed_file][feature.chrom] = feature.end - feature.start
    outfile.write(outheader)
    for chrom in cov_table[in_beds[0]]:
        outline = chrom
        for bed_file in in_beds:
            outline = '{}\t{}'.format(outline, cov_table[bed_file][chrom])
        outfile.write(outline)
    outfile.close()


@cli.command()
@click.option('--infoyaml', type=click.STRING,
              default='info.yaml',
              help='Yaml file which will be modified. Default: info.yaml')
@click.option('--force/--not-force',
              default=False,
              help='Forces addition of genomic information without checking '
                   'to see if index and fasta files exist on system. '
                   'Default: --not-force')
@click.option('--all/--main',
              default=False,
              help='Sets whether to include all chromosomes or just the main '
                   'ones. If a chromosome has "_" in the name, it will not be '
                   'included if main is set. Examples include '
                   'chromosomes with _random, _alt, and chrUn_ in the name. '
                   'Default: --main')
@click.argument('genome', type=click.STRING)
@click.argument('fasta', type=click.STRING)
@click.argument('index', type=click.STRING)
def add_genome(genome, fasta, index, infoyaml, force, all):
    """
    Adds genome information to info.yaml file.

    Takes in a genome name and appropriate information to info.yaml file.

    \b
    Required arguments:
    GENOME     UCSC name of genome
    INDEX      Path to BS Seeker2 index
    FASTA      Location of a fasta file containing all chromosomal sequences.
    """
    #Checks to see if necessary files exist on system
    if not force:
        assert os.path.exists(index), \
            'Failure: {} does not exist. \nPlease ensure you input the ' \
            'correct path or use the --force option.'.format(index)
        assert os.path.exists(fasta), \
            'Failure: {} does not exist. \nPlease ensure you input the ' \
            'correct path or use the --force option.'.format(fasta)

    workingdir = tempfile.mkdtemp()

    #Gets chromosome sizes by first fetching them using fetchChromSizes from
    #  UCSC genome browser web site. Then, parses the resulting file.
    chrom_sizes = os.path.join(workingdir, 'chroms.txt')
    fcs = resource_filename(wgbs_tools.__name__, '../external/fetchChromSizes')
    command = 'sh {} {} > {}'.format(fcs, genome, chrom_sizes)
    print('Running command: {}'.format(command))
    subprocess.check_call(command, shell=True)

    # Appends info.yaml file
    outfile = open(infoyaml, 'ab')
    outfile.write('{}:\n'.format(genome))
    outfile.write('  fasta: {}\n'.format(fasta))
    outfile.write('  index: {}\n'.format(index))
    outfile.write('  chroms:\n')
    with open(chrom_sizes, 'r') as cs_file:
        for line in cs_file:
            line_list = line.split('\t')
            if all:
                printline = '      {}: {}'.format(line_list[0], line_list[1])
                outfile.write(printline)
            else:
                if '_' not in line_list[0]:
                    printline = '      {}: {}'.format(line_list[0], line_list[1])
                    outfile.write(printline)
    outfile.close()
    shutil.rmtree(workingdir)
