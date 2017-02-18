"""
Utilities (or tools) used to modify, edit, or work with SAM format files.

These usually edit SAM files from BS Seeker 2. See link for more info:
https://github.com/BSSeeker/BSseeker2#output-format

"""

import subprocess
import logging
import re
from multiprocessing import Pool
import itertools
import pysam
from wgbs_tools import utilities
from multiprocessing import Process
import multiprocessing
import gzip


def bam_to_permeth(in_bam, out_prefix, bed_prefix, genome,
                   meth_type, strand_type, max_dup_reads, chroms, threads):
    """"""
    # Assert inputs are correct before starting everything
    assert meth_type in ['C', 'CG', 'CH', 'CHG', 'CHH'], \
        'ERROR! Methylation type needs to be C, CG, CH, CHG, or CHH. ' \
        'Methylation was set to: {}'.format(meth_type)
    assert strand_type in ['positive', 'negative' , 'both'], \
        'ERROR! Strand needs to be positive, negative, or both. ' \
        'Strand was set to: {}'.format(strand_type)
    # Multithread the processes
    pool = multiprocessing.Pool(threads)
    for chrom in chroms:
        chrom_length = chroms[chrom]
        out_bed = '{}_{}.bed.gz'.format(out_prefix, chrom)
        pool.apply_async(chr_bam_to_permeth,
                         args=(in_bam, out_bed, bed_prefix, genome, meth_type,
                               strand_type, max_dup_reads, chrom, chrom_length))
    pool.close()
    pool.join()


def chr_bam_to_permeth(in_bam, out_bed, bed_prefix, genome, meth_type,
                       strand_type, max_dup_reads, chrom, chrom_length):
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
    positions = {}
    count = 0

    # Previous read info
    prevstart = 0
    prevstrand = '+'
    dupcount = 1

    # Create compressed outfile and write header line
    outfile = gzip.open(out_bed, 'wb')
    header_line = 'track name={}{} description={}_{} useScore=0 itemRgb=On ' \
                  'db={}\n'.format(bed_prefix, chrom, bed_prefix, chrom, genome)
    outfile.write(header_line)

    # Open bam file and process each read, one at a time
    samfile = pysam.AlignmentFile(in_bam, 'rb')
    for read in samfile.fetch(chrom):
        start = read.reference_start
        methstring = read.get_tag('XM')
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
        if prevstart == start and prevstrand == strand:
            dupcount += 1
            if dupcount > max_dup_reads:
                continue
        else:
            dupcount = 1
            prevstart = start
            prevstrand = strand

        # Pulls out methylation information and adds it to methylation dict
        for search_char in search_chars:
            print(methstring, search_char)
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
