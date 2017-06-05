"""
Tests fastqtools module
"""

from wgbs_tools import fastqtools
import os
import subprocess


def test_qual_filter_fastq(qual_fastq, tmpdir, correct_qualfil_fastq):
    """Tests qual_filter_fastq function"""
    out_fastq = os.path.join(str(tmpdir), 'test_qualfilterout.fq')
    fastqtools.qual_filter_fastq(qual_fastq, out_fastq)
    # command = 'cp {} {}'.format(out_fastq, correct_qualfil_fastq)
    # subprocess.check_call(command, shell = True)
    with open(out_fastq, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_qualfil_fastq, \
            'Quality FASTQ filter error.'


def test_adapter_remove(adap_fastq, tmpdir,
                        correct_noadap_fastq, correct_adap_fastq):
    """Tests adapter_remove to ensure adapters gets removed properly and only
    desired reads remain"""
    noadap_fq = os.path.join(str(tmpdir), 'test_noadapout.fq')
    adaptrim_fq = os.path.join(str(tmpdir), 'test_adaptrim.fq')
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


def test_seq_to_searchseq1():
    """Tests seq_to_searchseq to convert all possible characters"""
    seq = 'ACGTRYSWKMBDHVN'
    correct_seq = 'ACGT[AG][CT][GC][AT][GT][AC][CGT][AGT][ACT][ACG][ACGT]'
    new_seq = fastqtools.seq_to_searchseq(seq)
    assert new_seq == correct_seq, 'Improper conversion of all possible ' \
                                   'characters to searchseq'


def test_seq_to_searchseq2():
    """Tests seq_to_searchseq to convert Line-1 characters"""
    seq = 'TTYGTGGTGYGTYGTTTTTTAAKTYG'
    correct_seq = 'TT[CT]GTGGTG[CT]GT[CT]GTTTTTTAA[GT]T[CT]G'
    new_seq = fastqtools.seq_to_searchseq(seq)
    assert new_seq == correct_seq, 'Improper conversion of Line1 to searchseq'


def test_meth_motif(line1_fastq, tmpdir):
    # Line-1 sequece (including modular K base not found in Line-1 pyroseq)
    seq = 'TTYGTGGTGYGTYGTTTTTTAAKTYG'
    out_fastq = os.path.join(str(tmpdir), 'test_line1.fq')
    correct_out = [0.55, 0.8, 0.4, 0.4, 0.6, 4, 5, 2, 5, 2, 5, 3, 5]
    results = fastqtools.meth_motif(line1_fastq, seq, out_fastq)
    assert results == correct_out, 'Methylation motif calculation not working'

