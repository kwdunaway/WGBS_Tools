"""
Utilities (or tools) used to modify, edit, or work with SAM format files.

These usually edit SAM files from BS Seeker 2. See link for more info:
https://github.com/BSSeeker/BSseeker2#output-format

"""

import subprocess
import logging
import re


def samsorted_to_permeth(in_sam, meth_type, out_prefix, bed_prefix, genome, 
                         meth_type, strand_type, max_dup_reads, chroms):
    #Load search_chars with chars based on meth_type
    if meth_type == 'CG':
        search_chars = ['X', 'x']
    elif meth_type == 'CH':
        search_chars = ['Z', 'z', 'Y', 'y']
    elif meth_type == 'CHG':
        search_chars = ['Y', 'y']
    elif meth_type == 'CHH':
        search_chars = ['Z', 'z']
    else:
        print('ERROR! meth_type = {}'.format(meth_type))
    if strand_type != 'positive' or strand_type != 'negative' or  \
            strand_type != 'combined':
        print('ERROR! strand_type = {}'.format(strand_type))

    # Initialize variables
    currentchrom = "Not Set Yet"
    methylation = {}
    positions = {}
    count = 0

    # Previous read info
    prevstart = 0
    prevstrand = "+"
    prevmethstring = ""
    dupcount = 1

    # Columns for formatting SAM files
    chrc = 2
    startc = 3
    methc = 14
    strandc = 11

    with open(in_sam, 'r') as sam_in_file:
        for line in sam_in_file:
            #TODO: Multiprocess this.

            #If header, skips
            if re.match(r'^@', line):
                continue

            #Splits line and defines variables
            #TODO: Confirm correct parsing of these variables
            line = line[:-1]
            line_list = line.split('\t')
            chrom = line_list[chrc]
            start = line_list[startc]
            methstring = line_list[methc][:5]
            strand = line_list[strandc][5:6]

            #Skips if read is not used with given parameters
            if chrom not in chroms:
                continue
            if strand == '-' and strand_type == 'positive':
                continue
            if strand == '+' and strand_type == 'negative':
                continue

            #If duplicate line, take information until max_dup_reads passed
            if prevstart == start and prevstrand == strand:
                dupcount += 1
                if dupcount > max_dup_reads:
                    continue
            else:
                dupcount = 1

            #If found search_chars, add to the methylation dict
            for searchchar in search_chars:
                if searchchar in prevmethstring:
                    #TODO: Create Addto_MethylationHash subroutine
                    # Addto_MethylationHash(\ % Methylation, $searchchars[0],
                        # $prevmethstring, $prevstart, $prevstrand)

            #On next chromosome, so print current one if not set yet
            if chrom != currentchrom:
                if currentchrom != 'Not Set Yet':
                    #TODO: Create Print_MethylationHash subroutine
                    #Print_MethylationHash(\ % methylation, $outprefix,
                        # $currentchrom, $bedprefix)
                # Reset variables
                methylation = {}
                currentchrom = chrom
                dupcount = 1
                logging.info('Starting {}'.format(chrom))

            # Increment
            prevmethstring = methstring
            prevstart = start
            prevstrand = strand

        # Finish last line/chromosome
        #Addto_MethylationHash(\ % Methylation, $searchchars[0],
            # $prevmethstring, $prevstart, $prevstrand)
        #Print_MethylationHash(\ % methylation, $outprefix,
            # $currentchrom, $bedprefix)




