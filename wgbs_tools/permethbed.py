"""
Functions that utilize Percentage Methylated Bed files (aka permeth.bed)

Utilizes pybedtools found at http://daler.github.io/pybedtools/

Note: If you use pybedtools in your work, please cite the pybedtools
manuscript and the BEDTools manuscript:
Dale RK, Pedersen BS, and Quinlan AR. 2011. Pybedtools: a flexible Python
library for manipulating genomic datasets and annotations. Bioinformatics
27(24):3423-3424.
Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for
comparing genomic features. Bioinformatics 26(6):841-842.
"""

import logging
import os
import gzip
from threading import Thread
from pybedtools import BedTool

from wgbs_tools import utilities

def meth_count(feature):
    """
    Returns the number of methylated reads from a permeth BedTools object.

    :param feature: feature line of a BedTools object line
    :return: int as the number of methylated reads
    """
    [perc, total] = utilities.show_value(feature.name).split('-')
    return int(float(perc)*float(total)+.5)


def total_count(feature):
    """
    Returns the number of reads from a permeth BedTools object line.

    :param feature: feature line of a BedTools object
    :return: int as the total number of reads covered
    """
    total = utilities.show_value(feature.name).split('-')[1]
    return int(total)


def chrom_meth(pm_sample, chrom, roi_chrom, mask, meth_dict):
    """
    Analyzes the methylation over the designated chromosome.

    :param pm_sample: string containing the sample name
    :param chrom: string containing the chromosome name
    :param roi_chrom: bedtools data structure containing the regions of
                      interest only over the designated chromosome.
    :param mask: bedtools data structure containing the masked areas of the
                 genome.
    :param meth_dict: dictionary directly written to containing all of the
                      ROI results for all samples over the designated chromsome.
    :return: Nothing (results are written directly to meth_dict)
    """
    permeth_name = '{}{}.bed'.format(pm_sample, chrom)
    if not os.path.exists(permeth_name):
        permeth_name = '{}{}.bed.gz'.format(pm_sample, chrom)
    logging.info('Processing %s.', extra=permeth_name)
    pm_full = BedTool(permeth_name)
    if mask[0].chrom == 'chrNONE':
        pm_masked = pm_full
    else:
        pm_masked = pm_full - mask
    pm_essential = pm_masked.intersect(roi_chrom, u=True)
    for roi_line in roi_chrom:
        start = int(roi_line.start)
        end = int(roi_line.end)
        meth = 0
        total = 0
        cpg = 0
        for pm_line in pm_essential.all_hits(roi_line):
            meth = meth + int(meth_count(pm_line))
            total = total + int(total_count(pm_line))
            cpg += 1
        meth_dict[start][end][pm_sample]['meth'] = int(meth)
        meth_dict[start][end][pm_sample]['total'] = int(total)
        meth_dict[start][end][pm_sample]['cpg'] = int(cpg)


def create_window_roi(window_roi, windowsize, chroms):
    """
    Creates a bed file with all of the window locations.

    :param window_roi: file name to output the windows file (a ROI bed file).
    :param windowsize: size of windows.
    :param chroms: dict containing chromosome names (keys) and lengths (values)
    :return: Nothing
    """
    windows = []
    for chrom in chroms:
        length = chroms[chrom]
        start = 1
        while start < length:
            end = start + windowsize - 1
            windows.append((chrom, start, end))
            start += windowsize
    bed_windows = BedTool(windows)
    bed_windows.saveas(window_roi)


def roi_meth(in_bed_prefixes, in_sample_list, out_table, mask_file, roi_file,
             min_read_count, min_cpg_count, min_file_count, raw_data_name,
             thread_count):
    """
    Creates a table with the methylation across desired Regions of
    Interest (ROI).

    :param in_bed_prefixes: list of bed file prefixes.
    :param in_sample_list: list of sample names. Order corresponds with
                           in_bed_prefixes.
    :param out_table: name of output table file.
    :param mask_file: bed or gtf file that will contain areas masked from
                      analysis (ie: any areas in this file will be ignored).
    :param roi_file: bed or gtf file containing the areas of the genome you
                     want analyzed.
    :param min_read_count: minimum read count necessary for a region of
                           interest. If a sample has less than this read count,
                           NA will be input instead of the average methylation
                           over the ROI.
    :param min_cpg_count: minimum CpG count necessary for a region of
                          interest. If a sample has less than this read count,
                          NA will be input instead of the average methylation
                          over the ROI.
    :param min_file_count: minimum file count to keep a region of interest. If
                           less than this many files/samples meet the
                           previous minimum requirements, that roi will not
                           have output in your out_table file.
    :param raw_data_name: optional file that (if populated) will be the output
                          of methylated and total read counts for each sample.
                          The minimums still apply and will work the same as
                          the main file.
    :param thread_count: int designating threads to allocate for multithreading.
    :return: Nothing
    """
    # Reduces thread count if there aren't enough tasks to fill all threads
    if len(in_bed_prefixes) < thread_count:
        thread_count = len(in_bed_prefixes)
    outfile = open(out_table, 'wb')
    header_line = 'chrom\tstart\tend\tname'
    for samp in in_sample_list:
        header_line = '{}\t{}'.format(header_line, samp)
    header_line = '{}\n'.format(header_line)
    outfile.write(header_line)
    print header_line
    if raw_data_name != "":
        raw_data = open(raw_data_name, 'wb')
        header_line = 'chrom\tstart\tend\tname'
        for samp in in_sample_list:
            header_line = '{0}\t{1}_methylated\t{1}_total\t{1}_cpgs'\
                .format(header_line, samp)
        header_line = '{}\n'.format(header_line)
        raw_data.write(header_line)

    roi = BedTool(roi_file)
    if mask_file != "":
        mask = BedTool(mask_file)
    else:
        mask = BedTool([('chrNONE', 0, 0)])

    # Get chromosome names in ROI file
    logging.info('Loading chromosomes:')
    chrom_names_tmp = []
    for line in roi:
        chrom = utilities.show_value(line.chrom)
        if chrom not in chrom_names_tmp:
            chrom_names_tmp.append(chrom)
    # Remove chromosome names without accompanying PerMeth file
    chrom_names = []
    for chrom in chrom_names_tmp:
        keepchrom = True
        for pm_sample in in_bed_prefixes:
            permeth_name = '{}{}.bed'.format(pm_sample, chrom)
            if not os.path.exists(permeth_name):
                permeth_name = '{}{}.bed.gz'.format(pm_sample, chrom)
                if not os.path.exists(permeth_name):
                    # logging.warning('Cannot access a file for {}, skipping!',
                    #                 extra=chrom)
                    print 'Cannot access a file {} for {}, skipping!'\
                        .format(permeth_name, chrom)
                    keepchrom = False
        if keepchrom:
            chrom_names.append(chrom)

    # Loop through, gather information, and print each chrom info
    for chrom in chrom_names:
        # Create methylation dictionary for chromosomal ROI
        roi_chrom = roi.all_hits(BedTool([(chrom, 0, 999999999)])[0])
        meth_dict = utilities.nested_dict(4, str)
        for feature in roi_chrom:
            meth_dict[feature.start][feature.end]['name'] = feature.name
        proc_list = list(in_bed_prefixes)
        def worker():
            """Worker for multithreading that analyzes a chromosome."""
            while proc_list:
                pm_prefix = proc_list.pop()
                chrom_meth(pm_prefix, chrom, roi_chrom, mask, meth_dict)
        threads = [Thread(target=worker) for i in range(thread_count)]
        [t.start() for t in threads]
        [t.join() for t in threads]

        # Print information into table
        for start in sorted(meth_dict):
            for end in sorted(meth_dict[start]):
                name = meth_dict[start][end]['name']
                print_line = '{}\t{}\t{}\t{}'.format(chrom, start, end, name)
                raw_col_line = print_line
                file_print_count = 0
                for pm_sample in in_bed_prefixes:
                    meth = meth_dict[start][end][pm_sample]['meth']
                    total = meth_dict[start][end][pm_sample]['total']
                    cpg = meth_dict[start][end][pm_sample]['cpg']
                    if total >= min_read_count and cpg >= min_cpg_count:
                        try:
                            float(meth)
                        except ValueError:
                            print "Not a float: {}".format(meth)
                        try:
                            float(total)
                        except ValueError:
                            print "Not a float: {}".format(total)
                        meth_perc = float(meth)/float(total)
                        print_line = '{0}\t{1:.3f}'.format(print_line, meth_perc)
                        file_print_count += 1
                    else:
                        print_line = '{0}\tNA'.format(print_line)
                    raw_col_line = '{}\t{}\t{}\t{}'\
                        .format(raw_col_line, meth, total, cpg)
                print_line = '{}\n'.format(print_line)
                raw_col_line = '{}\n'.format(raw_col_line)
                if file_print_count >= min_file_count:
                    outfile.write(print_line)
                    if raw_data_name != "":
                        raw_data.write(raw_col_line)


def convert_pm2dss(in_pmbed, out_dss):
    """
    Converts a single percent methylation bed file to DSS format.

    :param in_pmbed: percent methylation bed file name
    :param out_dss: DSS file name
    :return: Nothing
    """
    # Creates output DSS file and writes header line
    if out_dss.endswith('.gz'):
        dss = gzip.open(out_dss, 'wb')
    else:
        dss = open(out_dss, 'wb')
    print_line = 'chr\tpos\tN\tX\n'
    dss.write(print_line)

    # Open percent methylated bed file, process info, and prints to DSS file
    pm_essential = BedTool(in_pmbed)
    for pm_line in pm_essential:
        chrom = utilities.show_value(pm_line.chrom)
        start = int(pm_line.start)
        meth = int(meth_count(pm_line))
        total = int(total_count(pm_line))
        print_line = '{}\t{}\t{}\t{}\n'.format(chrom, start, total, meth)
        dss.write(print_line)
    dss.close()


def convert_pm2bg(in_pmbed, out_bg):
    """
    Converts a single percent methylation bed file to bedgraph format.

    :param in_pmbed: percent methylation bed file name
    :param out_bg: bedgraph file name
    :return: None
    """
    # Creates output bedgraph file and writes header line
    if out_bg.endswith('.gz'):
        bedgraph = gzip.open(out_bg, 'wb')
    else:
        bedgraph = open(out_bg, 'wb')

    # Open percent methylated bed file, process info, and prints to bg file
    pm_essential = BedTool(in_pmbed)
    for pm_line in pm_essential:
        chrom = utilities.show_value(pm_line.chrom)
        start = int(pm_line.start)
        end = start + 1
        perc = utilities.show_value(pm_line.name).split('-')[0]
        print_line = '{}\t{}\t{}\t{}\n'.format(chrom, start, end, perc)
        bedgraph.write(print_line)
    bedgraph.close()


def bed_meth_stats(in_bed):
    """
    Takes in a percent methylation bed file and returns the total CpG count
    as well as overall methylated read count and total read count. This is
    useful for determing conversion efficency as well as getting an overall
    methylation of a region or chromosome.

    :param in_bed: name of bed file
    :return: dictionary with methylation stats of the file
    """
    pm = BedTool(in_bed)
    meth = 0
    total = 0
    cpg = 0
    for pm_line in pm:
        meth = meth + int(meth_count(pm_line))
        total = total + int(total_count(pm_line))
        cpg += 1
    perc = float(meth) / float(total)
    meth_dict = {'perc': perc, 'meth': meth, 'total': total, 'cpgs': cpg}
    return meth_dict
