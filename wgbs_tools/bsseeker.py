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
    Aligns FASTQ file and outputs bam file using BS-Seeker2.

    :param bs2_path: Full path to bs_seeker2-align.py
    :param params: parameters passed to BS Seeker 2. See BS Seeker 2 manual for
    more details about all parameters.
    :param fasta: FASTA file of genome aligning reads to.
    :param bs2_index: BS Seeker 2 index of genome to align fastq to.
    :param fastq: Fastq file used as input for alignment.
    :param bam: BAM file to output results of BS-Seeker2
    :return:
    """
    #python BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3
    # -f bam -g /share/lasallelab/genomes/hg38/hg38.fa
    # -d /share/lasallelab/genomes/hg38/BSseek2_refgen/
    # -i raw_sequences/JLCM007A_noadap.fq.gz -o JLCM007A/JLCM007A_noadap.bam
    command = '{} {} -g {} -d {} -i {} -o {}'\
        .format(bs2_path, params, fasta, bs2_index, fastq, bam)
    logging.info(command)
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
    :return:
    """
    #python BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3
    # -f bam -g /share/lasallelab/genomes/hg38/hg38.fa
    # -d /share/lasallelab/genomes/hg38/BSseek2_refgen/
    # -i raw_sequences/JLCM007A_noadap.fq.gz -o JLCM007A/JLCM007A_noadap.bam
    command = '{} --aligner=bowtie2 {} -g {} -d {} -i {} -o {}'\
        .format(bs2_path, params, fasta, bs2_index, fq1, fq2, bam)
    logging.info(command)
    subprocess.check_call(command, shell=True)


def process_logs(noadap_log, trimmed_log, conv_eff_log, out_summary):
    """"""
    out_file = open(out_summary, 'wb')

    # Pull out alignment log information from log files
    usefulinfo = []
    for logfile in [noadap_log, trimmed_log]:
        with open(logfile, 'r') as conv_file:
            for linestr in conv_file:
                linestr = linestr.rstrip()
                line = linestr.split('\t')
                if bool(re.search('Read filename', linestr)):
                    usefulinfo.append(line[4])
                elif bool(re.search('Mappability', linestr)):
                    usefulinfo.append(line[3])
                elif bool(re.search('Total bases of uniquely mapped', linestr)):
                    usefulinfo.append(line[8])
                elif bool(re.search('mCG', linestr)):
                    usefulinfo.append(line[3])
                elif bool(re.search('mCHG', linestr)):
                    usefulinfo.append(line[3])
                elif bool(re.search('mCHH', linestr)):
                    usefulinfo.append(line[3])

    # Set up alignment log print line for sample information
    alignlogline =  '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'\
        .format(usefulinfo[0], usefulinfo[6], usefulinfo[2], usefulinfo[8],
                usefulinfo[1], usefulinfo[7], usefulinfo[3], usefulinfo[9],
                usefulinfo[4], usefulinfo[10], usefulinfo[5], usefulinfo[11])

    # Print out file with and without conversion efficency information
    if os.path.isfile(conv_eff_log):
        conv_file = open(conv_eff_log, 'r')
        cline = conv_file.next()
        cline = cline[:-1]
        print_line = 'FASTQ Files\t\tTotal bases of uniquely mapped reads' \
                 '\t\tMappability\t\tmCG\t\tmCHG\t\tmCHH\t{}'.format(cline)
        out_file.write(print_line)
        print_line = 'noadap\ttrimmed\tnoadap\ttrimmed\tnoadap\ttrimmed' \
                    '\tnoadap\ttrimmed\tnoadap\ttrimmed\tnoadap\ttrimmed'
        out_file.write(print_line)
        cline = conv_file.next()
        cline = cline[:-1]
        alignlogline = '{}\t{}'.format(alignlogline, cline)
        out_file.write(alignlogline)
        conv_file.close()
    else:
        print_line = 'FASTQ Files\t\tTotal bases of uniquely mapped reads' \
                 '\t\tMappability\t\tmCG\t\tmCHG\t\tmCHH'
        out_file.write(print_line)
        print_line = 'noadap\ttrimmed\tnoadap\ttrimmed\tnoadap\ttrimmed' \
                    '\tnoadap\ttrimmed\tnoadap\ttrimmed\tnoadap\ttrimmed'
        out_file.write(print_line)
        out_file.write(alignlogline)
    out_file.close()
