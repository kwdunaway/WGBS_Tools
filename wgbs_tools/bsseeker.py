"""
Commands relating to BS Seeker 2
see:
http://pellegrini.mcdb.ucla.edu/BS_Seeker2/

MIT License:
https://github.com/BSSeeker/BSseeker2/blob/master/LICENSE
"""

import subprocess
import logging

def align_bs2(bs2_path, params, fasta, bs2_index, fastq, bam):
    """

    :param bs2_path: Full path to bs_seeker2-align.py
    :param params: parameters passed to BS Seeker 2. See BS Seeker 2 manual for
    more details about all parameters.
    :param fasta: FASTA file of genome aligning reads
    :param bs2_index: BS Seeker 2 index
    :param fastq:
    :param bam:
    :return:
    """
    #python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3 -f bam -g /share/lasallelab/genomes/hg38/hg38.fa -d /share/lasallelab/genomes/hg38/BSseek2_refgen/ -i raw_sequences/JLCM007A_noadap.fq.gz -o JLCM007A/JLCM007A_noadap.bam
    command = '{} {} -g {} -d {} -i {} -o {}'\
        .format(bs2_path, params, fasta, bs2_index, fastq, bam)
    logging.info(command)
    subprocess.check_call(command, shell=True)



