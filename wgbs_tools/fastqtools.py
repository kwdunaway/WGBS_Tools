"""
Module containing all scripts that process and manipulate fastq files.
"""

import subprocess
import logging
import gzip
import multiprocessing
import re


def qual_filter_fastq(in_fastq, out_fastq):
    """
    Filters a fastq file based on the Y/N flag in the header line of each read.

    :param in_fastq: Input fastq file name. If ends in .gz, assumed gzipped.
                     Otherwise, it is assumed unzipped.
    :param out_fastq: Output fastq file name. If ends in .gz, it will zip the
                      resulting file. Otherwise, it is left unzipped.
    :return: String denoting completion of function or if the fastq file did
             not have a field for Illumina quality calls in the header lines.
    """
    # Checks to see if the fastq file has Illumina quality calls
    if in_fastq.endswith('.gz'):
        with gzip.open(in_fastq, 'r') as infq:
            first_line = infq.readline()
    else:
        with open(in_fastq, 'r') as infq:
            first_line = infq.readline()
    if ':Y:' not in first_line and ':N:' not in first_line:
        logging.info('no_qual_field')
        return

    # Filters fastq file
    if in_fastq.endswith(".gz"):
        if out_fastq.endswith(".gz"):
            command = 'gunzip -c {} | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' |   ' \
                      'grep -v "^--$" | gzip > {}'.format(in_fastq, out_fastq)
        else:
            command = 'gunzip -c {} | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' |   ' \
                      'grep -v "^--$" > {}'.format(in_fastq, out_fastq)
    else:
        if out_fastq.endswith(".gz"):
            command = 'cat {} | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' |   ' \
                      'grep -v "^--$" | gzip > {}'.format(in_fastq, out_fastq)
        else:
            command = 'cat {} | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' |   ' \
                      'grep -v "^--$" > {}'.format(in_fastq, out_fastq)
    #TODO: Enable logging
    # logging.info(command)
    print(command)
    subprocess.check_call(command, shell=True)


def adapter_remove(in_fastq, noadap_fq, adaptrim_fq, adap_seq, chew_length=10,
                   min_seqlength=35):
    """
    Separates fastq files for those with and without adapter sequence.

    If the reads did not contain adapter sequence, they are trimmed back
    chew_length bases. Otherwise, the reads that have adapter sequence are
    cut at the sequence and then trimmed back chew_length bases.

    :param in_fastq: input fastq file (can be .gz or uncompressed)
    :param noadap_fq: output fastq file name containing all reads without
                      adapter sequence in them. Trimmed back chew_length.
    :param adaptrim_fq: output fastq file name containing all reads that had
                        adapter sequence in them. The sequence was trimmed out
                        and trimmed back chew_length.
    :param adap_seq: Sequence of adapter that is used to filter for
                     adapter contamination.
    :param chew_length: Length of read that gets removed after adapter
                        sequence is found.
    :param min_seqlength: Minimum sequence length the read must be in order
                          to be kept. Any read shorter than this length will be
                          thrown out before writing to the adaptrim_fq file.
    :return: Nothing
    """
    # Opens write files (in either normal or compressed format)
    if noadap_fq.endswith('.gz'):
        noadap_outfile = gzip.open(noadap_fq, 'wb')
    else:
        noadap_outfile = open(noadap_fq, 'wb')
    if adaptrim_fq.endswith('.gz'):
        trimmed_outfile = gzip.open(adaptrim_fq, 'wb')
    else:
        trimmed_outfile = open(adaptrim_fq, 'wb')
    if in_fastq.endswith('.gz'):
        fastq_in_file = gzip.open(in_fastq, 'r')
    else:
        fastq_in_file = open(in_fastq, 'r')
    for header_line in fastq_in_file:
        seq_line = next(fastq_in_file)
        seq_line = seq_line[:-1]
        third_line = next(fastq_in_file)
        qual_line = next(fastq_in_file)
        qual_line = qual_line[:-1]
        seq_trimmed = seq_line.split(adap_seq)
        if len(seq_trimmed[0]) == len(seq_line):
            noadap_outfile.write(header_line)
            noadap_outfile.write('{}\n'.format(seq_line[:-chew_length]))
            noadap_outfile.write(third_line)
            noadap_outfile.write('{}\n'.format(qual_line[:-chew_length]))
        else:
            trimmed_seq_line = seq_trimmed[0][:-chew_length]
            if len(trimmed_seq_line) >= min_seqlength:
                qline_len = len(trimmed_seq_line)
                trimmed_outfile.write(header_line)
                trimmed_outfile.write('{}\n'.format(trimmed_seq_line))
                trimmed_outfile.write(third_line)
                trimmed_outfile.write('{}\n'.format(qual_line[:qline_len]))
    fastq_in_file.close()
    noadap_outfile.close()
    trimmed_outfile.close()


def pe_adapter_remove(fin_fastq, fnoadap_fq, fadaptrim_fq, fadap_seq,
                      rin_fastq, rnoadap_fq, radaptrim_fq, radap_seq,
                      chew_length, min_seqlength, threads):
    """
    Separates fastq files for those with and without adapter sequence.

    If the reads did not contain adapter sequence, they are trimmed back
    chew_length bases. Otherwise, the reads that have adapter sequence are
    cut at the sequence and then trimmed back chew_length bases.

    :param fin_fastq: forward input fastq file (can be .gz or uncompressed)
    :param fnoadap_fq: forward output fastq file name containing all reads
                       without adapter sequence in them. Trimmed back
                       chew_length.
    :param fadaptrim_fq: forward output fastq file name containing all reads
                         that had adapter sequence in them. The sequence was
                         trimmed out and trimmed back chew_length.
    :param fadap_seq: Sequence of forward adapter that is used to filter for
                      adapter contamination.
    :param rin_fastq: reverse input fastq file (can be .gz or uncompressed)
    :param rnoadap_fq: reverse output fastq file name containing all reads
                       without adapter sequence in them. Trimmed back
                       chew_length.
    :param radaptrim_fq: reverse output fastq file name containing all reads
                         that had adapter sequence in them. The sequence was
                         trimmed out and trimmed back chew_length.
    :param radap_seq: Sequence of reverse adapter that is used to filter for
                      adapter contamination.
    :param chew_length: Length of read that gets removed after adapter
                        sequence is found.
    :param min_seqlength: Minimum sequence length the read must be in order
                          to be kept. Any read shorter than this length will be
                          thrown out before writing to the adaptrim_fq file.
    :return: Nothing
    """
    # Opens write files (in either normal or compressed format)
    if fnoadap_fq.endswith('.gz'):
        fnoadap_outfile = gzip.open(fnoadap_fq, 'wb')
    else:
        fnoadap_outfile = open(fnoadap_fq, 'wb')
    if fadaptrim_fq.endswith('.gz'):
        ftrimmed_outfile = gzip.open(fadaptrim_fq, 'wb')
    else:
        ftrimmed_outfile = open(fadaptrim_fq, 'wb')
    if rnoadap_fq.endswith('.gz'):
        rnoadap_outfile = gzip.open(rnoadap_fq, 'wb')
    else:
        rnoadap_outfile = open(rnoadap_fq, 'wb')
    if radaptrim_fq.endswith('.gz'):
        rtrimmed_outfile = gzip.open(radaptrim_fq, 'wb')
    else:
        rtrimmed_outfile = open(radaptrim_fq, 'wb')

    if fin_fastq.endswith('.gz'):
        for_file = gzip.open(fin_fastq, 'r')
    else:
        for_file= open(fin_fastq, 'r')
    if rin_fastq.endswith('.gz'):
        rev_file = gzip.open(rin_fastq, 'r')
    else:
        rev_file= open(rin_fastq, 'r')

    for for_header in for_file:
        # Get forward read information
        for_seq = next(for_file)
        for_seq = for_seq[:-1]
        for_third = next(for_file)
        for_qual = next(for_file)
        for_qual = for_qual[:-1]
        for_seq_trimmed = for_seq.split(fadap_seq)

        # Get reverse read information
        rev_header = next(rev_file)
        rev_seq = next(rev_file)
        rev_seq = rev_seq[:-1]
        rev_third = next(rev_file)
        rev_qual = next(rev_file)
        rev_qual = rev_qual[:-1]
        rev_seq_trimmed = rev_seq.split(radap_seq)

        cutlength = len(for_seq_trimmed[0])
        if len(rev_seq_trimmed[0]) < cutlength:
            cutlength = len(rev_seq_trimmed[0])

        pattern = re.compile('^@.* [^:]*:Y:[^:]*:')
        if pattern.match(for_header) or pattern.match(rev_header):
            print('Bad read: {}'.format(for_header))
        elif cutlength == len(for_seq):
            fnoadap_outfile.write(for_header)
            fnoadap_outfile.write('{}\n'.format(for_seq[:-chew_length]))
            fnoadap_outfile.write(for_third)
            fnoadap_outfile.write('{}\n'.format(for_qual[:-chew_length]))
            rnoadap_outfile.write(rev_header)
            rnoadap_outfile.write('{}\n'.format(rev_seq[:-chew_length]))
            rnoadap_outfile.write(rev_third)
            rnoadap_outfile.write('{}\n'.format(rev_qual[:-chew_length]))
        else:
            cutlength -= chew_length
            if cutlength >= min_seqlength:
                for_trimmed_seq = for_seq[:cutlength]
                ftrimmed_outfile.write(for_header)
                ftrimmed_outfile.write('{}\n'.format(for_trimmed_seq))
                ftrimmed_outfile.write(for_third)
                ftrimmed_outfile.write('{}\n'.format(for_qual[:cutlength]))
                rev_trimmed_seq = rev_seq[:cutlength]
                rtrimmed_outfile.write(rev_header)
                rtrimmed_outfile.write('{}\n'.format(rev_trimmed_seq))
                rtrimmed_outfile.write(rev_third)
                rtrimmed_outfile.write('{}\n'.format(rev_qual[:cutlength]))
    for_file.close()
    rev_file.close()
    fnoadap_outfile.close()
    ftrimmed_outfile.close()
    rnoadap_outfile.close()
    rtrimmed_outfile.close()
