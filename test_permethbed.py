"""
Tests permethbed module(s)
"""

from wgbs_tools import permethbed
import os
import gzip
import sys


def test_meth_count(feature):
    """Tests meth_count"""
    meth = permethbed.meth_count(feature)
    assert meth == 5, 'meth_count is not returning methylation count correctly'


def test_total_count(feature):
    """Tests meth_count"""
    meth = permethbed.meth_count(feature)
    assert meth == 5, 'meth_count is not returning methylation count correctly'


def test_roi_meth(bed_folder, working_dir):
    """Tests roi_meth to multithread analysis"""
    in_bed_prefixes = []
    for num in range(6):
        fileprefix = os.path.join(bed_folder, 'pm0{}_'.format(num))
        in_bed_prefixes.append(fileprefix)
    in_sample_list = ['pm01', 'pm02', 'pm03', 'pm04', 'pm05', 'pm06']
    roi_file = os.path.join(bed_folder, 'roi1.bed')
    mask_file = ''
    # out_table = os.path.join(working_dir, 'roi_meth_table1.txt')
    out_table = 'roi_meth_table1.txt'
    out_2col_name = ""
    min_read_count = 0
    min_file_count = 0
    thread_count = 1
    permethbed.roi_meth(in_bed_prefixes, in_sample_list, out_table, mask_file,
                        roi_file, min_read_count, min_file_count,
                        out_2col_name, thread_count)