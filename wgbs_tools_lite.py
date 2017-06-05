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
from wgbs_tools import samutils
from wgbs_tools import permethbed
from wgbs_tools import fastqtools

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
                   'Default: 1')
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
@click.option('--verbose', default=False, is_flag=True)
@click.argument('input_tsv', type=click.STRING)
@click.argument('out_table', type=click.STRING)
@click.argument('roi_file', type=click.STRING)
def roi(input_tsv, out_table, roi_file, mask, min_read_count, min_cpg_count,
        min_sample_coverage, raw_data, threads, verbose):
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
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
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
    permethbed.roi_meth(in_bed_prefixes, in_sample_list, out_table,
                        mask, roi_file, min_read_count, min_cpg_count,
                        min_sample_coverage, raw_data, threads)


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
                   'Default: 1')
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
              default=default_info_yaml,
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. '
                   'Default: info.yaml')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('input_tsv', type=click.STRING)
@click.argument('out_table', type=click.STRING)
def window(input_tsv, out_table, windowsize, mask, raw_data, threads,
           min_read_count, min_cpg_count, min_sample_coverage, genome,
           infoyaml, verbose):
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
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
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
              default=False,
              help='Boolean which indicates if the output files will be '
                   'compressed (.gz format) or not. Default: --no-gz')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def pm2dss(in_prefix, out_prefix, gz, verbose):
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
@click.option('--suffix', type=click.STRING,
              default='.bed.gz',
              help='Requires all in files end in a particular string. '
                   'Default: .bed.gz')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_file', type=click.STRING)
def pm_stats(in_prefix, out_file, suffix, verbose):
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
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    out = open(out_file, 'wb')
    outline = 'name\tpercentage\tmethylated_reads\ttotal_reads\tcpg_count\n'
    out.write(outline)
    for in_file_name in glob.glob('{}*{}'.format(in_prefix, suffix)):
        uniqname = in_file_name.split(in_prefix)[1].split(suffix)[0]
        meth_dict = permethbed.bed_meth_stats(in_file_name)
        outline = '{}\t{}\t{}\t{}\t{}\n' \
            .format(uniqname, meth_dict['perc'], meth_dict['meth'],
                    meth_dict['total'], meth_dict['cpgs'])
        out.write(outline)
    out.close()


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
                   'Default: 1')
@click.option('--header', 'header_name', type=click.STRING,
              default='',
              help='Name in header of bed file. This is useful for loading '
                   'these files into a genome browser. The name of the track '
                   'will be {header}_chr#. '
                   'Default: bed_prefix of lowest directory')
@click.option('--infoyaml', type=click.STRING,
              default=default_info_yaml,
              help='Yaml file which contains information which could change '
                   'based on experiment. Read README.md to modify the '
                   'default or create your own. '
                   'Default: info.yaml')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_bam', type=click.STRING)
@click.argument('bed_prefix', type=click.STRING)
def bam2pm(in_bam, bed_prefix, genome, methtype, strand, max_dup_reads, threads,
           header_name, infoyaml, verbose):
    """
    Converts BAM to pm_bed.

    Converts a BAM file (produced by either BS_Seeker2 or Bismark) into a
    Percentage Methylation BED format (pm_bed).

    \b
    Required arguments:
    IN_BAM       Input sorted BAM file.
                 It should be indexed as well, although this will index the
                 bam file if an index is not found.
    BED_PREFIX   Prefix of output percent methylation bed files.
                 If you want to output in a directory other than the current
                 working directory, use Option: out_dir.
    """
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

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

    # Convert sam to permeth bed files (percent methylation bed files)
    # This also only prints CpGs on chromosomes in ranges defined in the yaml
    samutils.bam_to_permeth(in_bam, bed_prefix, header_prefix, genome,
                            methtype, strand, max_dup_reads, chroms, threads)


@cli.command()
@click.option('--infoyaml', type=click.STRING,
              default=default_info_yaml,
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
@click.option('--verbose', default=False, is_flag=True)
@click.argument('genome', type=click.STRING)
@click.argument('fasta', type=click.STRING)
@click.argument('index', type=click.STRING)
def add_genome(genome, fasta, index, infoyaml, force, all, verbose):
    """
    Adds genome information to info.yaml file.

    Takes in a genome name and appropriate information to info.yaml file.

    \b
    Required arguments:
    GENOME     UCSC name of genome
    INDEX      Path to BS Seeker2 index
    FASTA      Location of a fasta file containing all chromosomal sequences.
    """
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    #Checks to see if necessary files exist on system
    if not force:
        assert os.path.exists(index), \
            'Failure: {} does not exist. \nPlease ensure you input the ' \
            'correct_output path or use the --force option.'.format(index)
        assert os.path.exists(fasta), \
            'Failure: {} does not exist. \nPlease ensure you input the ' \
            'correct_output path or use the --force option.'.format(fasta)

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


@cli.command()
@click.option('--verbose', default=False, is_flag=True)
@click.argument('in_prefix', type=click.STRING)
@click.argument('out_prefix', type=click.STRING)
def pm2bg(in_prefix, out_prefix, gz, verbose):
    """
    Converts pm_bed to bedgraph format.

    The percent methylation bed files (pm_bed) are used to create a bedGraph
    file, which can be used for various reasons, including creating a trackhub.

    \b
    Required arguments:
    IN_PREFIX      Prefix for all percent methlyated bed files (can be
                   compressed and/or uncompressed)
    OUT_prefix     Prefix for all outfiles (in bedgraph format). The full
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
    # TODO: convert this to take in multiple permeth files and output a single bg
    if gz:
        suffix = '.bg.gz'
    else:
        suffix = '.bg'
    for bed_name in glob.glob('{}*.bed'.format(in_prefix)):
        uniqname = bed_name.split(in_prefix)[1].split('.bed')[0]
        bg_name = '{}{}{}'.format(out_prefix, uniqname, suffix)
        logging.info('Processing {}'.format(bed_name))
        logging.info('Creating {}'.format(bg_name))
        permethbed.convert_pm2bg(bed_name, bg_name)
    for bed_name in glob.glob('{}*.bed.gz'.format(in_prefix)):
        uniqname = bed_name.split(in_prefix)[1].split('.bed.gz')[0]
        bg_name = '{}{}{}'.format(out_prefix, uniqname, suffix)
        logging.info('Processing {}'.format(bed_name))
        logging.info('Creating {}'.format(bg_name))
        permethbed.convert_pm2bg(bed_name, bg_name)


@cli.command()
@click.option('--seq', type=click.STRING,
              default='',
              help='')
@click.option('--table/--not-table',
              default=False,
              help='The input file as a 2 column table of fastq files instead '
                   'of a single fastq file name. '
                   'Default: --not-table')
@click.option('--working-dir', 'working_dir', type=click.STRING,
              default='',
              help='Working directory where temp files are written and read. '
                   'Default: <EMPTY> (Uses tempfile.mkdtemp() to define a '
                   'temperary working directory)')
@click.option('--verbose', default=False, is_flag=True)
@click.argument('fastq_file', type=click.STRING)
@click.argument('out_table', type=click.STRING)
@click.argument('seq', type=click.STRING)
def motif(in_file, out_table, seq, list, working_dir, verbose):
    """
    Adds genome information to info.yaml file.

    Takes in a genome name and appropriate information to info.yaml file.

    \b
    Required arguments:
    IN_FILE           Input fastq file.
                      If --table is used, the file is expected to be a 2
                      column table of fastq files instead. Ex:
                          sample1.fq.gz  sample1
                          sample2.fq     sample2
    OUT_TABLE         Name of output table.
    SEQ               Sequence motif to search for. For example:
                      Human Line1: TTYGTGGTGYGTYGTTTTTTAAKTYGGTT
    """
    if working_dir == '':
        workingdir = tempfile.mkdtemp()
    else:
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        workingdir = working_dir
    if verbose:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    #set global variables
    printcols = 0
    num_cpgs = seq.count('Y')

    #create out table and header within it
    out_file = open(out_table, 'wb')
    out_file.write('Name\tAvg_Total')
    for n in range(1, num_cpgs+1):
        out_file.write('\tAvg_CpG_{}'.format(n))
        printcols += 1
    for n in range(1, num_cpgs+1):
        out_file.write('\tCpG_meth_{}'.format(n))
        out_file.write('\tCpG_total_{}'.format(n))
        printcols += 2
    out_file.write('\n')

    #Checks to see if the list flag is true or false
    if list:
        #If there is a list of fastq files
        list_file = open(in_file, 'r')
        for line in list_file:
            line = line.rstrip()
            fastq_file = line.split('\t')[0]
            fastq_name = line.split('\t')[1]
            out_fastq = os.path.join(workingdir, '{}.fq'.format(fastq_name))
            meth = fastqtools.meth_motif(fastq_file, seq, out_fastq)

            #print out methylation in table
            out_file.write('{}'.format(fastq_name))
            for num in range(0, printcols):
                out_file.write('\t{}'.format(meth[num]))
            out_file.write('\n')

    else:
        # If there is just one fastq file
        fastq_file = in_file
        fastq_name = os.path.split(fastq_file)[1]
        out_fastq = os.path.join(workingdir, '{}.fq'.format(fastq_name))
        meth = fastqtools.meth_motif(fastq_file, seq, out_fastq)

        # print out methylation in table
        out_file.write('{}'.format(fastq_name))
        for num in range(0, printcols):
            out_file.write('\t{}'.format(meth[num]))
        out_file.write('\n')

    out_file.close()

    #Remove temp files
    if working_dir == '':
        shutil.rmtree(workingdir)


