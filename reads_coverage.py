#!/usr/bin/env python3
# coding=utf-8
'''
Program to statistics of read coverage.
Usage: python reads_coverage.py bamfile out_prefix out_bam_file
bamfile: nanopore DRS fasta/fastq aligned to transcriptome with minimap2(minimap2 -ax map-ont). 
out_bam_file: output file(Bam file containing all full-length reads) path.
Tips: The read that spanning at least 95% of the length their best hit of an existing transcript will be defined full-length reads.

Author: lt
Date: 2021-01-09 15:26:47
LastEditTime: 2021-04-02 09:55:54
'''

import pysam
import sys,csv,os

try:
    bamfile = sys.argv[1]
    out_prefix = sys.argv[2]
    out_bam_file = sys.argv[3]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

def path_prefix(out_prefix):
    """
    Generate output path based on prefix.
    """
    prefix = os.path.basename(out_prefix)
    stat = prefix + '.stat.txt'
    data = prefix + '.data.csv'

    if os.path.dirname(out_prefix) == '':
        stat_out_path = os.path.join(os.getcwd(),stat)
        data_out_path = os.path.join(os.getcwd(),data)
    else:
        stat_out_path = os.path.join(os.path.dirname(out_prefix),stat)
        data_out_path = os.path.join(os.path.dirname(out_prefix),data)
    
    return stat_out_path,data_out_path


def check_read(read):
    """
    callback function for total reads
    """
    if read.is_secondary or read.is_supplementary:
        return False        
    else:
        return True

def coverage(bamfile,out_prefix,out_bam_file):
    """Statistics of read coverage(only calculate the primary alignment), will output 2 file.
    Args:
        bamfile : Transcriptome-aligned sequence
        out_prefix : Output file name prefix
    
    Returns:
        prefix.stat : Percentage and number of full-length reads
        prefix.data : Coverage and reference name of each read
    """
    full_length_reads = 0
    header = ['reference_name','reads_name','reads_coverage']
    out_path = path_prefix(out_prefix)

    with open(out_path[0],'w',newline='') as stat,open(out_path[1],'w',newline='') as data:
        read_data = csv.writer(data,delimiter='\t')
        read_data.writerow(header)
        bamfile = pysam.AlignmentFile(bamfile,'rb')
        full_length_bam_file = pysam.AlignmentFile(out_bam_file, "wb", template=bamfile) #creat full_length bam file.
        read_number = bamfile.count(until_eof=True,read_callback=check_read)
        for reads in bamfile.fetch():
            if reads.is_supplementary or reads.is_secondary:
                continue
            else:
                reference_name = reads.reference_name
                reads_name = reads.query_name
                reads_coverage = reads.reference_length / bamfile.get_reference_length(reference_name)
                read_data.writerow([reference_name,reads_name,reads_coverage])
                if reads_coverage >= 0.95:
                    full_length_reads += 1
                    full_length_bam_file.write(reads)
        bamfile.close()
        full_length_bam_file.close()
        percentage = full_length_reads / read_number
        stat.write('Number of reads representing full-length transcripts:{}{}'.format(full_length_reads,'\n'))
        stat.write('Percentage of reads representing full-length transcripts:{}'.format(percentage))

if __name__ == '__main__':
    
    coverage(bamfile,out_prefix,out_bam_file) 
    
