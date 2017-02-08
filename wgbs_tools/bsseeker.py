"""
Functions that wrap BS-Seeker2 calls.

MIT License of BS-Seeker2:
https://github.com/BSSeeker/BSseeker2/blob/master/LICENSE

For more information on BS-Seeker2, see documenation at:
http://pellegrini.mcdb.ucla.edu/BS_Seeker2/

"""

import subprocess
import logging


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
