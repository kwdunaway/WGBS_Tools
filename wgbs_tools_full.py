"""
Click wrapper to run all of the wgbs_tools much easier and user friendly.
"""

import logging
logger = logging.getLogger(__name__)
import os
import glob
import subprocess
import tempfile
import shutil
import click
import yaml
import wgbs_tools

from pkg_resources import resource_filename
from wgbs_tools import fastqtools
from wgbs_tools import bsseeker
from wgbs_tools import samutils
from wgbs_tools import permethbed
import wgbs_tools_lite as lite


# Allows the default yaml file to be in WGBS_Tools directory
default_info_yaml = resource_filename(wgbs_tools.__name__, '../info.yaml')

# Default for number of CPUs is set to 1 because most people want to test
# code as default
NUM_CPUS = 1
# If you want to have the CPU default be the number of CPUs available,
# uncomment the following lines:
#import multiprocessing
#NUM_CPUS = multiprocessing.cpu_count()


@click.group()
def cli():
    """A versatile toolkit to manipulate and analyze WGBS data."""
    pass


# Adds all of the commands from the lite version
cli.add_command(lite.add_genome)
cli.add_command(lite.bam2pm)
cli.add_command(lite.pm2dss)
cli.add_command(lite.pm_stats)
cli.add_command(lite.roi)
cli.add_command(lite.window)


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
                   'Default: 1')
@click.option('--working-dir', 'working_dir', type=click.STRING,
              default='',
              help='Working directory where temp files are written and read. '
                   'Default: <EMPTY> (Uses tempfile.mkdtemp() to define a '
                   'temperary working directory)')
@click.option('--infoyaml', type=click.STRING,
              default=default_info_yaml,
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. '
                   'Default: info.yaml')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_fastq', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def process_se(in_fastq, out_prefix, out_dir, chew_length, min_seqlength,
               genome, noadap_bs2_params, adaptrim_bs2_params, methtype,
               strand, max_dup_reads, conv_chrom, threads, working_dir,
               infoyaml, verbose):
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
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    if working_dir == '':
        workingdir = tempfile.mkdtemp()
    else:
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        workingdir = working_dir
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
    samutils.bam_to_permeth(full_bam, bed_prefix, out_prefix, genome,
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
                    conv_eff_dict['cpgs'])
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
                   'Default: 1')
@click.option('--working_dir', type=click.STRING,
              default='',
              help='Working directory where temp files are written and read. '
                   'Default: <EMPTY> (Uses tempfile.mkdtemp() to define a '
                   'temperary working directory)')
@click.option('--infoyaml', type=click.STRING,
              default=default_info_yaml,
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. '
                   'Default: info.yaml')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_fastq_f', type=click.STRING)
@click.argument('in_fastq_r', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def process_pe(in_fastq_f, in_fastq_r, out_prefix, out_dir, genome, chew_length,
               min_readlength, noadap_bs2_params, adaptrim_bs2_params,
               methtype, strand, max_dup_reads, conv_chrom, threads,
               working_dir, infoyaml, verbose):
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
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
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
    samutils.bam_to_permeth(full_bam, bed_prefix, out_prefix, genome,
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
@click.option('--verbose', default=False, is_flag=True)
@click.argument('noadap_log', type=click.STRING)
@click.argument('adaptrim_log', type=click.STRING)
@click.argument('out_summary', type=click.STRING)
def sumlogs(noadap_log, adaptrim_log, out_summary, conv_eff, verbose):
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
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    bsseeker.process_logs(noadap_log, adaptrim_log, conv_eff, out_summary)

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
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def adjustcols(in_prefix, out_prefix, suffix, cols, adjusts, header, verbose):
    """
    Adjusts numerical column of files.

    Takes all files matching the prefix of in_prefix and outputs

    \b
    Required arguments:
    IN_PREFIX    Prefix of all files input into the
    OUT_PREFIX   Prefix of all output files
    """
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
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
                   'Default: 1')
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
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_fastq', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def trim_sefq(in_fastq, out_prefix, adapter, out_dir, threads, chew_length,
              min_seqlength, working_dir, verbose):
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
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    if working_dir == '':
        workingdir = tempfile.mkdtemp()
    else:
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        workingdir = working_dir
    if out_dir != '':
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

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
                   'Default: 1')
@click.option('--chew', 'chew_length', type=click.INT,
              default=10,
              help='Number of bases to removed off of the end of each read. '
                   'Default: 10')
@click.option('--min-read', 'min_seqlength', type=click.INT,
              default=35,
              help='Minimum read length of reads after adapter trimming and '
                   'chew. If the read is less than this lenght, it is not '
                   'included in the output. Default: 35')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_for_fq', type=click.STRING)
@click.argument('in_rev_fq', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def trim_pefq(in_for_fq, in_rev_fq, out_prefix, for_adap, rev_adap, out_dir,
              threads, chew_length, min_seqlength, verbose):
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
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

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
