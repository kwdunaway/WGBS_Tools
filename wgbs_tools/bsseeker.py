"""
Functions that wrap BS-Seeker2 calls.

MIT License of BS-Seeker2:
https://github.com/BSSeeker/BSseeker2/blob/master/LICENSE

For more information on BS-Seeker2, see documenation at:
http://pellegrini.mcdb.ucla.edu/BS_Seeker2/
"""

import subprocess
import logging
import os
import re


def align_bs2(bs2_path, params, fasta, bs2_index, fastq, bam):
    """
    Aligns single end FASTQ file and outputs bam file using BS-Seeker2.

    This is a simple wrapper that allows customization of BS-Seeker2 calling.
    The options given allow the user to call extra parameters if they wish, but
    requires essential values to be used.

    An example BS-Seeker2 call:

    python BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3
     -f bam -g /share/lasallelab/genomes/hg38/hg38.fa
     -d /share/lasallelab/genomes/hg38/BSseek2_refgen/
     -i raw_sequences/JLCM007A_noadap.fq.gz -o JLCM007A/JLCM007A_noadap.bam

    :param bs2_path: Full path to bs_seeker2-align.py
    :param params: parameters passed to BS Seeker 2. See BS Seeker 2 manual for
    more details about all parameters.
    :param fasta: FASTA file of genome aligning reads to.
    :param bs2_index: BS Seeker 2 index of genome to align fastq to.
    :param fastq: Fastq file used as input for alignment.
    :param bam: BAM file to output results of BS-Seeker2
    :return: Nothing
    """
    command = '{} {} -g {} -d {} -i {} -o {}'\
        .format(bs2_path, params, fasta, bs2_index, fastq, bam)
    #TODO: Enable logging
    # logging.info(command)
    print(command)
    subprocess.check_call(command, shell=True)


def align_bs2_pe(bs2_path, params, fasta, bs2_index, fq1, fq2, bam):
    """
    Aligns FASTQ file and outputs bam file using BS-Seeker2.

    :param bs2_path: Path to bs_seeker2-align.py (can either be full path or
                     call from $PATH)
    :param params: parameters passed to BS Seeker 2. See BS Seeker 2 manual for
                   more details about all parameters.
    :param fasta: FASTA file of genome aligning reads to.
    :param bs2_index: BS Seeker 2 index of genome to align fastq to.
    :param fastq: Fastq file used as input for alignment.
    :param bam: BAM file to output results of BS-Seeker2
    :return: Nothing
    """
    #python BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3
    # -f bam -g /share/lasallelab/genomes/hg38/hg38.fa
    # -d /share/lasallelab/genomes/hg38/BSseek2_refgen/
    # -i raw_sequences/JLCM007A_noadap.fq.gz -o JLCM007A/JLCM007A_noadap.bam
    command = '{} --aligner=bowtie2 {} -g {} -d {} -i {} -o {}'\
        .format(bs2_path, params, fasta, bs2_index, fq1, fq2, bam)
    #TODO: Enable logging
    # logging.info(command)
    print(command)
    subprocess.check_call(command, shell=True)


def process_logs(noadap_log, trimmed_log, conv_eff_log, out_summary):
    """
    Processes basic information from BS-Seeker2 logs.

    This results in a summary file of the following information from the
    BS-Seeker2 log files:
      1) Percentage of reads that were able to map to the genome
      2) Number of bases
      3) CG methylation percentage
      4) CHG methylation percentage
      5) CHH methylation percentage

    :param noadap_log: log file for the no_adapter contamination alignment
    :param trimmed_log: log file for the adapter trimmed alignment
    :param conv_eff_log: log file for conversion efficiency. It is possible
                         for this file not to exist. If that is the case,
                         this function will modify the output summary file to
                         not include those columns.
    :param out_summary: output summary file
    :return: Nothing
    """
    out_file = open(out_summary, 'wb')

    # Pull out alignment log information from log files
    usefulinfo = []
    for logfile in [noadap_log, trimmed_log]:
        with open(logfile, 'r') as conv_file:
            for linestr in conv_file:
                linestr = linestr.rstrip()
                line = linestr.split()
                if bool(re.search('Read filename', linestr)):
                    usefulinfo.append(line[4])
                elif bool(re.search('Mappability', linestr)):
                    usefulinfo.append(line[4])
                elif bool(re.search('Total bases of uniquely mapped', linestr)):
                    usefulinfo.append(line[9])
                elif bool(re.search('mCG', linestr)):
                    usefulinfo.append(line[3])
                elif bool(re.search('mCHG', linestr)):
                    usefulinfo.append(line[3])
                elif bool(re.search('mCHH', linestr)):
                    usefulinfo.append(line[3])


    # Sets up the print line with desired information in correct orientation
    alignlogline =  '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'\
        .format(usefulinfo[2], usefulinfo[8],
                usefulinfo[1], usefulinfo[7], usefulinfo[3], usefulinfo[9],
                usefulinfo[4], usefulinfo[10], usefulinfo[5], usefulinfo[11])

    # Print out file with and without conversion efficency information
    first_line = 'Total bases of uniquely mapped reads' \
                 '\t\tMappability\t\tmCG\t\tmCHG\t\tmCHH'
    second_line = 'noadap\ttrimmed\tnoadap\ttrimmed' \
                  '\tnoadap\ttrimmed\tnoadap\ttrimmed\tnoadap\ttrimmed\n'

    # If the conversion efficiency file exists
    if os.path.isfile(conv_eff_log):
        conv_file = open(conv_eff_log, 'r')
        cline = conv_file.next()
        first_line = '{}\t{}'.format(first_line, cline)
        out_file.write(first_line)
        out_file.write(second_line)
        cline = conv_file.next()
        cline = cline[:-1]
        alignlogline = '{}\t{}\n'.format(alignlogline, cline)
        out_file.write(alignlogline)
        conv_file.close()

    # If the conversion efficiency file does not exist
    else:
        first_line = '{}\n'.format(first_line)
        out_file.write(first_line)
        out_file.write(second_line)
        alignlogline = '{}\n'.format(alignlogline)
        out_file.write(alignlogline)

    out_file.close()
