"""
Defines test variables
"""

import wgbs_tools
import pytest
import tempfile
import os
import gzip
import shutil
from pkg_resources import resource_filename
import yaml
from pybedtools import BedTool
import sys

# working_directory = tempfile.mkdtemp()

# @pytest.yield_fixture(scope='session')
# def temp_dir():
#     tmpdir = tempfile.mkdtemp()
#     def closure(*names):
#         path = os.path.join(tmpdir, *names)
#         try:
#             os.mknod(path)
#         except OSError as exc:
#             if exc.errno != 17:
#                 raise
#         return path
#     yield closure
#     shutil.rmtree(tmpdir)
#

# @pytest.fixture
# def working_dir():
#     """Get a working directory for all files"""
#     return working_directory

@pytest.fixture
def sys_version():
    """Returns system major number for test_utilities"""
    return sys.version_info.major

@pytest.fixture
def bs2_index():
    """
    Path of test bs2 index for test_align_bs2.

    This is a fake genome that is small and useful only for testing to see if
    the pipeline works. Useful only for testing code purposes, not for other
    errors derived outside the scope of this package.
    """
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/genome')

@pytest.fixture
def fasta():
    """Fasta file for fake genome for test_align_bs2"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/genome.fa')

@pytest.fixture
def chroms():
    """Chromosomes of fake genome as a dict with key=name and value=size for
    test_align_bs2"""
    chroms = {'chr1': 1000, 'chr2': 500}

@pytest.fixture
def bs2_path():
    """Gets information from info.yaml in root directory of package and returns
    the path to BS Seeker 2 for test_align_bs2"""
    infoyaml = resource_filename(wgbs_tools.__name__, '../info.yaml')
    stream = file(infoyaml, 'r')
    info_dict = yaml.safe_load(stream)
    return info_dict['bs2_path']

@pytest.fixture
def bs2_index():
    """BSSeeker2 index for small genome"""
    return resource_filename(wgbs_tools.__name__, '../tests/data/small_genome/')

@pytest.fixture
def noadap_fastq():
    """FASTQ file for test_align_bs2"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/noadap.fq')

@pytest.fixture
def correct_noadapsam():
    """Correct output for test_align_bs2_noadap in test_align_bs2"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/noadap.sam')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def trimmed_fastq():
    """FASTQ file for test_align_bs2"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/trimmed.fq')

@pytest.fixture
def correct_trimmedsam():
    """Correct output for test_align_bs2_adaptrim in test_align_bs2"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/trimmed.sam')
    # return path
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def qual_fastq():
    """Input FASTQ for test_qual_filter_fastq"""
    return resource_filename(wgbs_tools.__name__, '../tests/data/fastq/qual.fq')

@pytest.fixture
def correct_qualfil_fastq():
    """Correct outpur for test_qual_filter_fastq"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/fastq/correctfil.fq')
    # return path
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def adap_fastq():
    """Input FASTQ for test_adapter_remove"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/fastq/adapterreads.fq')

@pytest.fixture
def correct_noadap_fastq():
    """Correct no adapter contamination file output for test_adapter_remove"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/fastq/correctnoadap.fq')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_adap_fastq():
    """Correct adapter contamination trim file output for test_adapter_remove"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/fastq/correctadap.fq')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def test_bam():
    """Bam file for testing samutils"""
    return resource_filename(wgbs_tools.__name__, '../tests/data/bam/test.bam')

@pytest.fixture
def testpe_bam():
    """Bam file for testing samutils"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/bam/test_pe.bam')

@pytest.fixture
def correct_chr1bed():
    """Correct chr1 percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bam/correct_chr1.bed')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_chr2bed():
    """Correct chr2 percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bam/correct_chr2.bed')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_chrNHbed():
    """Correct chrNH percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bam/correct_chrNH.bed')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correctpe_chr1bed():
    """Correct chr1 percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bam/correctpe_1.bed')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correctpe_chr2bed():
    """Correct chr2 percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bam/correctpe_2.bed')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correctpe_chrNHbed():
    """Correct chrNH percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bam/correctpe_chrNH.bed')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def feature():
    """Represents a single line of a bed file. Used in test_permethbed"""
    return BedTool([('chrF', 0, 1, '0.5-10')])[0]

@pytest.fixture
def bed_folder():
    """Folder containing all bedfiles"""
    path = resource_filename(wgbs_tools.__name__, '../tests/data/bed/')
    return path

@pytest.fixture
def correct_meth_table1():
    """Correct chrF percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/roi_meth_table1.txt')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_meth_raw1():
    """Correct chrF percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/roi_meth_table1.raw.txt')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_meth_table2():
    """Correct chrF percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/roi_meth_table2.txt')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_meth_raw2():
    """Correct chrF percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/roi_meth_table2.raw.txt')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_window_roi():
    """Correct content for making window ROI"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/window_roi.bed')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_meth_table3():
    """Correct chrF percent methylation bed file content for tests in
    test_samtutils"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/roi_meth_table3.txt')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_out_dss():
    """Correct output for DSS conversion"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/pm01_chrF.dss')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_out_bg():
    """Correct output for bedgraph conversion"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/bed/pm01_chrF.bg')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def line1_fastq():
    """FASTQ containing Line1 reads to be analyzed with test_meth_motif"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/fastq/line1.fq')
    return path
