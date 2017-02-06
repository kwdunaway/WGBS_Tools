"""
Tests fastqtools module
"""

from wgbs_tools import fastqtools
import os
import subprocess


def test_qual_filter_fastq(qual_fastq, working_dir, correct_qualfil_fastq):
    """Tests qual_filter_fastq function"""
    out_fastq = os.path.join(working_dir, 'test_qualfilterout.fq')
    fastqtools.qual_filter_fastq(qual_fastq, out_fastq)
    # command = 'cp {} {}'.format(out_fastq, correct_qualfil_fastq)
    # subprocess.check_call(command, shell = True)
    with open(out_fastq, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_qualfil_fastq, \
            'Quality FASTQ filter error.'


def test_adapter_remove(adap_fastq, working_dir,
                        correct_noadap_fastq, correct_adap_fastq):
    """Tests adapter_remove to ensure adapters gets removed properly and only
    desired reads remain"""
    noadap_fq = os.path.join(working_dir, 'test_noadapout.fq')
    adaptrim_fq = os.path.join(working_dir, 'test_adaptrim.fq')
    adap_seq = 'AGATCGGAAG'
    chew_length = 10
    min_seqlength = 35
    fastqtools.adapter_remove(adap_fastq, noadap_fq, adaptrim_fq, adap_seq,
                              chew_length, min_seqlength)
    with open(noadap_fq, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_noadap_fastq, \
            'Adapter removal error in printing no adapter contamination output.'
    with open(adaptrim_fq, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_adap_fastq, \
            'Adapter removal error in printing adapter contamination output.'


