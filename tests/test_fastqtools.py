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


