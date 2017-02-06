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

working_directory = tempfile.mkdtemp()

@pytest.yield_fixture(scope='session')
def temp_dir():
    tmpdir = tempfile.mkdtemp()
    def closure(*names):
        path = os.path.join(tmpdir, *names)
        try:
            os.mknod(path)
        except OSError as exc:
            if exc.errno != 17:
                raise
        return path
    yield closure
    shutil.rmtree(tmpdir)


@pytest.fixture
def working_dir():
    """Get a working directory for all files"""
    return working_directory

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
    return resource_filename(wgbs_tools.__name__, '../tests/data/qual.fq')

@pytest.fixture
def correct_qualfil_fastq():
    """Correct outpur for test_qual_filter_fastq"""
    path = resource_filename(wgbs_tools.__name__, '../tests/data/correctfil.fq')
    # return path
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def adap_fastq():
    """Input FASTQ for test_adapter_remove"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/adapterreads.fq')

@pytest.fixture
def correct_noadap_fastq():
    """Correct no adapter contamination file output for test_adapter_remove"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/correctnoadap.fq')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

@pytest.fixture
def correct_adap_fastq():
    """Correct adapter contamination trim file output for test_adapter_remove"""
    path = resource_filename(wgbs_tools.__name__,
                             '../tests/data/correctadap.fq')
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content
