"""
Utilities (or tools) used to modify, edit, or work with SAM format files.

These usually edit SAM files from BS Seeker 2. See link for more info:
https://github.com/BSSeeker/BSseeker2#output-format

"""

import multiprocessing
import gzip
import pysam
from wgbs_tools import utilities

# import subprocess
# import logging
# import re
# from multiprocessing import Pool
# import itertools
# from multiprocessing import Process


def bam_to_permeth(in_bam, out_prefix, bed_prefix, genome,
                   meth_type, strand_type, max_dup_reads, chroms, threads):
    """

    :param in_bam: input bam file name
    :param out_prefix: output prefix for bed files
    :param bed_prefix: header prefix within bed files produced
    :param genome: genome name (ex: mm10 or hg38)
    :param meth_type: methylation type (ex: CG or CH)
    :param strand_type: string that indicates strandedness. If you are
                        interested in strand specific methyation, change this to
                        '+' or '-'. Otherwise, use 'both'
    :param max_dup_reads: int of maximum duplicate reads. For lower coverage
                          WGBS experiments, use 1. If there are more than
                          this amount of reads with the same start and end
                          locations, the extras will be thrown out (and
                          assumed to be PCR duplicates).
    :param chroms: dict of chromosome names (keys) and lengths (values)
    :param threads: int of number of threads for multiprocessing
    :return: Nothing
    """
    # Assert inputs are correct before starting everything
    assert meth_type in ['C', 'CG', 'CH', 'CHG', 'CHH'], \
        'ERROR! Methylation type needs to be C, CG, CH, CHG, or CHH. ' \
        'Methylation was set to: {}'.format(meth_type)
    assert strand_type in ['positive', 'negative', 'both'], \
        'ERROR! Strand needs to be positive, negative, or both. ' \
        'Strand was set to: {}'.format(strand_type)
    # Multithread the processes
    #TODO: Allow both BS_seeker2 and Bismark inputs. Currently only works for
    #  Bismark
    pool = multiprocessing.Pool(threads)
    for chrom in chroms:
        chrom_length = chroms[chrom]
        out_bed = '{}{}.bed.gz'.format(out_prefix, chrom)
        pool.apply_async(chr_bam_to_permeth,
                         args=(in_bam, out_bed, bed_prefix, genome, meth_type,
                               strand_type, max_dup_reads, chrom, chrom_length))
    pool.close()
    pool.join()


def chr_bam_to_permeth(in_bam, out_bed, bed_prefix, genome, meth_type,
                       strand_type, max_dup_reads, chrom, chrom_length):
    """

    :param in_bam:
    :param out_bed:
    :param bed_prefix:
    :param genome:
    :param meth_type:
    :param strand_type:
    :param max_dup_reads:
    :param chrom:
    :param chrom_length:
    :return:
    """
    #TODO: Incorporate chrom_length
    if meth_type == 'C':
        search_chars = ['Z', 'z', 'Y', 'y', 'X', 'x']
    elif meth_type == 'CG':
        search_chars = ['X', 'x']
    elif meth_type == 'CH':
        search_chars = ['Z', 'z', 'Y', 'y']
    elif meth_type == 'CHG':
        search_chars = ['Y', 'y']
    elif meth_type == 'CHH':
        search_chars = ['Z', 'z']
    else:
        raise ValueError('ERROR! Methylation type needs to be C, CG, CH, CHG, '
                         'or CHH. Methylation was set to: {}'.format(meth_type))
    if strand_type == 'positive':
        analyzed_strands = ['+']
    elif strand_type == 'negative':
        analyzed_strands = ['-']
    elif strand_type == 'both':
        analyzed_strands = ['+', '-']
    else:
        raise ValueError('ERROR! Strand needs to be positive, negative, or '
                         'both. Strand was set to: {}'.format(strand_type))

    # Initialize variables
    methylation = {}
    # positions = {}
    # count = 0

    # Previous read info
    prevstart = 0
    prevstrand = '+'
    dupcount = 1
    prevmpos = 0

    # Create compressed outfile and write header line
    outfile = gzip.open(out_bed, 'wb')
    header_line = 'track name={}{} description={}{} useScore=0 itemRgb=On ' \
                  'db={}\n'.format(bed_prefix, chrom, bed_prefix, chrom, genome)
    outfile.write(header_line)

    # Open bam file and process each read, one at a time
    samfile = pysam.AlignmentFile(in_bam, 'rb')
    for read in samfile.fetch(str(chrom)):
        start = read.reference_start
        methstring = read.get_tag('XM')
        mpos = read.mpos
        if mpos > read.reference_start:
            methlen = mpos - read.reference_start
            methstring = methstring[0:methlen]
        strand = read.get_tag('XO')[0]
        if strand == '+':
            strandmult = 1
        elif strand == '-':
            strandmult = -1
        else:
            raise ValueError('ERROR! Strand from bam file was {}, not the '
                             'expected + or -'.format(strand))

        # Skips if read is not on correct strand
        if strand not in analyzed_strands:
            continue
        # If duplicate line, take information until max_dup_reads passed
        if prevstart == start and prevstrand == strand and prevmpos == mpos:
            dupcount += 1
            if dupcount > max_dup_reads:
                continue
        else:
            dupcount = 1
            prevstart = start
            prevstrand = strand
            prevmpos = mpos

        # Pulls out methylation information and adds it to methylation dict
        for search_char in search_chars:
            offsets = utilities.find_occurences(methstring, search_char)
            for pos in offsets:
                basepos = start + pos * strandmult
                if basepos in methylation:
                    methylation[basepos] = \
                        '{}{}'.format(methylation[basepos], search_char)
                else:
                    methylation[basepos] = search_char
    for start in sorted(methylation.iterkeys()):
        methstring = methylation[start]
        end = start + 1
        meth = sum(1 for c in methstring if c.isupper())
        total = len(methstring)
        meth_perc = float(meth)/float(total)
        meth_field = '{0:.2f}-{1}'.format(meth_perc, total)
        if meth_perc == 0:
            color = '0,0,0'  # black (if meth = 0)
        elif meth_perc <= .6:
            color = '27,74,210'  # blue (if 0 < meth <= .6)
        elif meth_perc <= .8:
            color = '27,210,57'  # green (if .6 < meth <= .8)
        else:
            color = '210,27,27'  # red (if meth > .8)
        bed_line = '{}\t{}\t{}\t{}\t0\t+\t0\t0\t{}\n'\
            .format(chrom, start, end, meth_field, color)
        outfile.write(bed_line)
    outfile.close()
