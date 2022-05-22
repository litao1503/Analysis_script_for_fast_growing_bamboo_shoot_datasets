#!/usr/bin/env python3
# coding=utf-8
'''
Program to statistics of read length.

Usage: python reads_statistics.py input_bam lname > out.stat
input_bam: accept bam file,and if there are multiple bams they can be provided in a comma-separated list.
lname: sample name comma-separated list.
Important: recommend use no filter bam file.

Author: lt
Date: 2021-01-08 10:03:16
LastEditTime: 2021-01-19 10:57:57
'''

import pysam
import numpy as np
import pandas as pd
import sys 

try:
    input_bam = sys.argv[1]
    lname = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

def N50(lengths):
    
    """
    Abstract: Returns the N50 length of the list of reads length. 
    Usage:    N50(numlist)
    Based on the dinovski/n50.py:
    https://gist.github.com/dinovski/2bcdcc770d5388c6fcc8a656e5dbe53c
    """
    #sort reads_length longest>shortest
    reads_length = sorted(lengths,reverse=True)
    csum = np.cumsum(reads_length)
    median = sum(lengths)/2
    
    # get index for cumsum >= N/2
    csum_median = min(csum[csum >= median])
    ind = np.where(csum == csum_median)
    
    n50 = reads_length[int(ind[0])]  
    return n50    

def STAT(read_length,mapped_length,mapped_read_length,read_number,match_length):

    average_length = np.mean(read_length)
    median_length = np.median(read_length)
    N50_length = N50(read_length)
    maximun_read_length = max(read_length)
    average_mapped_length = np.mean(mapped_length) 
    median_mapped_length = np.median(mapped_length)
    maximun_mapped_length = max(mapped_length)
    N50_mapped_length = N50(mapped_length)
    average_reads_accuracy = sum(match_length) / sum(mapped_read_length)
    reads_mapped_ratio = sum(mapped_length) / sum(mapped_read_length)
    
    
    stat = [read_number,average_length,median_length,N50_length,maximun_read_length,average_mapped_length,
            median_mapped_length,maximun_mapped_length,N50_mapped_length,average_reads_accuracy,reads_mapped_ratio]
    
    return stat   

def ReadStatistics(bam_file):
    
    read_number = 0
    read_length = []
    mapped_length = []
    mapped_read_length = []
    match_length = []

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    for reads in bamfile.fetch(until_eof=True):
        if reads.is_supplementary or reads.is_secondary:
            continue
        else:
            read_number += 1
            read_length.append(reads.query_length)
            if reads.is_unmapped == False:
                mapped_length.append(reads.query_alignment_length)
                mapped_read_length.append(reads.query_length)
                match_length.append(reads.get_cigar_stats()[0][0])  # get all match base(CIGAR = M) number

    
    stat = STAT(read_length,mapped_length,mapped_read_length,read_number,match_length)

    return stat,read_length,mapped_length,mapped_read_length,read_number,match_length
    


if __name__ == '__main__':
    #read number is the row read number, which is equal read number in alignment input fasta file.
    index = ['Read number','Average length','Median length','N50 length','Maximun read length','Average mapped length',
             'Median mapped length','Maximun mapped length','N50 mapped length','Average reads accuracy','Reads mapped ratio']
    samples = []
    reads_stat = pd.DataFrame(index=index)
    all_reads_length = []
    all_mapped_length = []
    all_mapped_read_length = []
    all_read_number = 0
    all_match_length = []
    for bam,name in zip(input_bam.split(','),lname.split(',')): 
        readstatistics = ReadStatistics(bam)
        reads_stat[name] = readstatistics[0]
        samples.append(name)
        all_reads_length += readstatistics[1]
        all_mapped_length += readstatistics[2]
        all_mapped_read_length += readstatistics[3]
        all_read_number += readstatistics[4]
        all_match_length += readstatistics[5]
    reads_stat['Overall'] = STAT(all_reads_length,all_mapped_length,all_mapped_read_length,all_read_number,all_match_length) 
    reads_stat.to_csv(sys.stdout,sep='\t')
