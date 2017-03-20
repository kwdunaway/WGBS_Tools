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
    total = permethbed.total_count(feature)
    assert total == 10, 'total_count is not returning total count correctly'


def test_create_window_roi(tmpdir, correct_window_roi):
    """Tests the create_window_roi function"""
    window_roi = os.path.join(str(tmpdir), 'window_roi.bed')
    # window_roi = 'window_roi.bed'
    windowsize = 10
    chroms = {'chrF': 200}
    permethbed.create_window_roi(window_roi, windowsize, chroms)
    with open(window_roi, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_window_roi, \
            'Error in creating correct_output window roi.'


def test_roi_meth1(bed_folder, tmpdir, correct_meth_table1,
                  correct_meth_raw1):
    """Tests roi_meth to multithread analysis"""
    in_bed_prefixes = []
    for num in range(6):
        numname = num + 1
        fileprefix = os.path.join(bed_folder, 'pm0{}_'.format(numname))
        in_bed_prefixes.append(fileprefix)
    in_sample_list = ['pm01', 'pm02', 'pm03', 'pm04', 'pm05', 'pm06']
    roi_file = os.path.join(bed_folder, 'roi1.bed')
    mask_file = ''
    out_table = os.path.join(str(tmpdir), 'roi_meth_table1.txt')
    out_raw_name = os.path.join(str(tmpdir), 'roi_meth_table1.raw.txt')
    min_read_count = 1
    min_cpg_count = 1
    min_file_count = 1
    thread_count = 2
    permethbed.roi_meth(in_bed_prefixes, in_sample_list, out_table, mask_file,
                        roi_file, min_read_count, min_cpg_count, min_file_count,
                        out_raw_name, thread_count)
    with open(out_table, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_meth_table1, \
            'Error in creating correct_output methylation table from roi.'
    with open(out_raw_name, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_meth_raw1, \
            'Error in creating correct_output raw methylation table from roi.'


def test_roi_meth2(bed_folder, tmpdir, correct_meth_table2,
                  correct_meth_raw2):
    """Tests roi_meth to mask analysis correctly"""
    in_bed_prefixes = []
    for num in range(6):
        numname = num + 1
        fileprefix = os.path.join(bed_folder, 'pm0{}_'.format(numname))
        in_bed_prefixes.append(fileprefix)
    in_sample_list = ['pm01', 'pm02', 'pm03', 'pm04', 'pm05', 'pm06']
    roi_file = os.path.join(bed_folder, 'roi1.bed')
    mask_file = 'tests/data/bed/mask2.bed'
    out_table = os.path.join(str(tmpdir), 'roi_meth_table2.txt')
    out_raw_name = os.path.join(str(tmpdir), 'roi_meth_table2.raw.txt')
    min_read_count = 1
    min_cpg_count = 1
    min_file_count = 1
    thread_count = 1
    permethbed.roi_meth(in_bed_prefixes, in_sample_list, out_table, mask_file,
                        roi_file, min_read_count, min_cpg_count, min_file_count,
                        out_raw_name, thread_count)
    with open(out_table, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_meth_table2, \
            'Error in creating correct_output methylation table from roi.'
    with open(out_raw_name, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_meth_raw2, \
            'Error in creating correct_output raw methylation table from roi.'


def test_roi_meth3(bed_folder, tmpdir, correct_meth_table3):
    """Tests roi_meth to mask analysis correctly across windows"""
    in_bed_prefixes = []
    for num in range(6):
        numname = num + 1
        fileprefix = os.path.join(bed_folder, 'pm0{}_'.format(numname))
        in_bed_prefixes.append(fileprefix)
    in_sample_list = ['pm01', 'pm02', 'pm03', 'pm04', 'pm05', 'pm06']
    roi_file = os.path.join(bed_folder, 'window_roi.bed')
    mask_file = 'tests/data/bed/mask2.bed'
    out_table = os.path.join(str(tmpdir), 'roi_meth_table3.txt')
    out_raw_name = ''
    min_read_count = 1
    min_cpg_count = 1
    min_file_count = 1
    thread_count = 1
    permethbed.roi_meth(in_bed_prefixes, in_sample_list, out_table, mask_file,
                        roi_file, min_read_count, min_cpg_count, min_file_count,
                        out_raw_name, thread_count)
    with open(out_table, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_meth_table3, \
            'Error in creating correct_output methylation table from roi.'


def test_convert_pm2dss(bed_folder, tmpdir, correct_out_dss):
    """Tests the converter of pm to dss formats"""
    in_pmbed = os.path.join(bed_folder, 'pm01_chrF.bed')
    out_dss = os.path.join(str(tmpdir), 'pm01_chrF.dss')
    permethbed.convert_pm2dss(in_pmbed, out_dss)
    with open(out_dss, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_out_dss, \
            'Error in converting permeth to DSS format.'


def test_convert_pm2bg(bed_folder, tmpdir, correct_out_bg):
    """Tests the converter of pm to bg formats"""
    in_pmbed = os.path.join(bed_folder, 'pm01_chrF.bed')
    out_dss = os.path.join(str(tmpdir), 'pm01_chrF.bg')
    permethbed.convert_pm2bg(in_pmbed, out_dss)
    with open(out_dss, 'r') as content_file:
        testcontent = content_file.read()
        assert testcontent == correct_out_bg, \
            'Error in converting permeth to DSS format.'


def test_bed_meth_stats(bed_folder):
    """Tests bed_meth_stats"""
    in_pmbed = os.path.join(bed_folder, 'pm01_chrF.bed')
    results_dict = permethbed.bed_meth_stats(in_pmbed)
    correct_dict = {'perc': 0.5, 'meth': 130, 'total': 260, 'cpgs': 5}
    assert results_dict == correct_dict, 'Error in getting bed meth stats.'
