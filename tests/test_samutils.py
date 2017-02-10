"""
Tests samutils module(s)
"""

from wgbs_tools import samutils
import os
import gzip


def test_chr_bam_to_permeth(test_bam, working_dir, correct_chr1bed):
    """Tests bam_to_permeth in samutils"""
    chrom = 'chr1'
    bed_prefix = 'test'
    out_bed = '{}{}_{}.bed.gz'.format(working_dir, bed_prefix, chrom)
    genome = 'mm10'
    meth_type = 'CG'
    strand_type = 'combined'
    max_dup_reads = 1
    chrom_length = 99999999
    samutils.chr_bam_to_permeth(test_bam, out_bed, bed_prefix, genome,
                                meth_type, strand_type, max_dup_reads, chrom,
                                chrom_length)
    with gzip.open(out_bed, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_chr1bed, \
            'BAM to Percent Methylation conversion is not working correctly.'


def test_bam_to_permeth(test_bam, working_dir, correct_chr1bed,
                        correct_chr2bed, correct_chrNHbed):
    """Tests multiprocessing of sam to permethbed conversion"""
    chroms = {'chr1': 99999999, 'chr2': 99999999, 'chrNH': 5}
    bed_prefix = 'test'
    out_prefix = os.path.join(working_dir, bed_prefix)
    genome = 'mm10'
    meth_type = 'CG'
    strand_type = 'combined'
    max_dup_reads = 1
    threads = 1
    samutils.bam_to_permeth(test_bam, out_prefix, bed_prefix, genome,
                   meth_type, strand_type, max_dup_reads, chroms, threads)
    out_bed = '{}_chr1.bed.gz'.format(out_prefix)
    with gzip.open(out_bed, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_chr1bed, \
            'BAM to Percent Methylation conversion is not working correctly.'
    out_bed = '{}_chr2.bed.gz'.format(out_prefix)
    with gzip.open(out_bed, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_chr2bed, \
            'BAM to Percent Methylation conversion is not working correctly ' \
            'for multiprocessing.'
    out_bed = '{}_chrNH.bed.gz'.format(out_prefix)
    with gzip.open(out_bed, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_chrNHbed, \
            'Multiprocess permeth creation not working for empty chromosome ' \
            'information.'
