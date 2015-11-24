# This is a model bash script for running analysis on 
# whole genome bisulfite sequencing data.             
# Parameter names and function calls need to be 
# changed to fit each user's needs.
# Absolute file pathing is recommended when possible.

# Set file pathing
PATH=$PATH:/share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/

# Load toolkits such as Bowtie, SAMtools, pysam, etc.
module load bowtie/1.1.1
module load samtools
module load sratoolkit
module load bedtools2
export PYTHONPATH=/share/lasallelab/pysam/lib/python2.7/site-packages/

# Example Run in genome: mm10 for sample name: JLAC001A

gunzip -c raw_sequences/Sample_JLAC001A/*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > raw_sequences/JLAC001A_filtered.fq
gzip raw_sequences/JLAC001A_filtered.fq

# Create two sets of the data with no adapter and with trimmed
perl /share/lasallelab/programs/perl_script/adapter_split.pl raw_sequences/JLAC001A_filtered.fq.gz raw_sequences/JLAC001A_noadap.fq.gz raw_sequences/JLAC001A_withadap.fq.gz
perl /share/lasallelab/programs/perl_script/adapter_trimmer.pl raw_sequences/JLAC001A_withadap.fq.gz raw_sequences/JLAC001A_trimmed.fq.gz 45 10

# Align sequences to mm10
mkdir JLAC001A
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 23 -e 90 -m 3 -f bam -g /share/lasallelab/genomes/mm10/mm10.fa -d /share/lasallelab/genomes/mm10/BSseek2_refgen/ -i raw_sequences/JLAC001A_noadap.fq.gz -o JLAC001A/JLAC001A_noadap.bam --unmapped=JLAC001A/JLAC001A_noadap_unmapped.fa --multiple-hit=JLAC001A/JLAC001A_noadap_multihit.fa
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 23 -e 80 -m 2 -f bam -g /share/lasallelab/genomes/mm10/mm10.fa -d /share/lasallelab/genomes/mm10/BSseek2_refgen/ -i raw_sequences/JLAC001A_trimmed.fq.gz -o JLAC001A/JLAC001A_trimmed.bam --unmapped=JLAC001A/JLAC001A_trimmed_unmapped.fa --multiple-hit=JLAC001A/JLAC001A_trimmed_multihit.fa

# Sort and convert results from bam files to sam files
samtools sort JLAC001A/JLAC001A_noadap.bam JLAC001A/JLAC001A_noadap_sorted
samtools sort JLAC001A/JLAC001A_trimmed.bam JLAC001A/JLAC001A_trimmed_sorted
samtools merge JLAC001A/JLAC001A.bam JLAC001A/JLAC001A_noadap_sorted.bam JLAC001A/JLAC001A_trimmed_sorted.bam
samtools view JLAC001A/JLAC001A.bam > JLAC001A/JLAC001A.sam

# Take sam output and convert to bed files with percentage methylation
mkdir JLAC001A/tmp
perl /share/lasallelab/programs/perl_script/SAMsorted_to_permeth.pl JLAC001A/JLAC001A.sam JLAC001A/tmp/PerMeth_JLAC001A JLAC001A mm10 CG combined

# Ensure that the results can be viewed on the genome browser and are valid
mkdir JLAC001A/PerMeth_JLAC001A
perl /share/lasallelab/programs/perl_script/gbcompliance.pl mm10 JLAC001A/tmp/PerMeth_JLAC001A_ JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_ JLAC001A JLAC001A
rm -r JLAC001A/tmp

# Remove the positions also found in CpG islands
mkdir JLAC001A/NoCGI_Permeth_JLAC001A
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr1.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr1.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr2.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr2.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr3.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr3.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr4.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr4.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr5.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr5.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr6.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr6.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr7.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr7.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr8.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr8.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr9.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr9.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr10.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr10.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr11.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr11.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr12.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr12.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr13.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr13.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr14.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr14.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr15.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr15.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr16.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr16.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr17.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr17.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr18.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr18.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr19.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr19.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr20.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr20.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr21.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr21.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chr22.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chr22.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chrX.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chrX.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chrY.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chrY.bed
bedtools subtract -a JLAC001A/PerMeth_JLAC001A/PerMeth_JLAC001A_chrM.bed -b /share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed > JLAC001A/NoCGI_Permeth_JLAC001A/NoCGI_Permeth_JLAC001A_chrM.bed


