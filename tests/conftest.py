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


# @pytest.fixture
# def working_dir():
#     """Get a working directory for all files"""
#     return working_directory


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

def fasta():
    """Fasta file for fake genome for test_align_bs2"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/genome.fa')

def chroms():
    """Chromosomes of fake genome as a dict with key=name and value=size for
    test_align_bs2"""
    chroms = {'chr1': 1000, 'chr2': 500}

def bs2_path():
    """Gets information from info.yaml in root directory of package and returns
    the path to BS Seeker 2 for test_align_bs2"""
    infoyaml = resource_filename(wgbs_tools.__name__, '../info.yaml')
    stream = file(infoyaml, 'r')
    info_dict = yaml.safe_load(stream)
    return info_dict['bs2_path']

def bs2_fastq():
    """FASTQ file for test_align_bs2"""
    return resource_filename(wgbs_tools.__name__,
                             '../tests/data/small_genome/test.fq')

def noadap_bam():
    """Output of bam file created through test_align_bs2"""
    return os.path.join(working_directory, 'noadap_test_out.bam')

def adaptrim_bam():
    """Output of bam file created through test_align_bs2"""
    return os.path.join(working_directory, 'adaptrim_test_out.bam')


