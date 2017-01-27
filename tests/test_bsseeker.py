"""
Tests bsseeker
"""

from wgbs_tools import bsseeker
import os
import subprocess


def test_align_bs2_noadap(bs2_path, fasta, bs2_index, bs2_fastq, noadap_bam):
    """
    Tests function align_bs2 from bsseeker using reads without adapter
    contamination trimmed out.
    """
    params = '-m 3 -f bam'
    bsseeker.align_bs2(bs2_path, params, fasta, bs2_index, bs2_fastq,
                       noadap_bam)


def test_align_bs2_adaptrim(bs2_path, fasta, bs2_index, bs2_fastq,
                            adaptrim_bam):
    """
    Tests function align_bs2 from bsseeker using reads with adapter
    contamination trimmed out of the reads, as well as 10bp chewed back.
    """
    params = '-m 2 -f bam'
    bsseeker.align_bs2(bs2_path, params, fasta, bs2_index, bs2_fastq,
                       adaptrim_bam)

