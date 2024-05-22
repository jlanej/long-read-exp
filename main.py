# utilities for long read sequencing data analysis
import os
import sys
import pysam
import argparse

utility_list = ['create_aligned_fasta']

# extract all the reads from the alignment file to a list
def extract_reads(alignment_file):
    a=pysam.AlignmentFile(alignment_file, "rb")
    reads = []
    for read in a.fetch():
        reads.append(read)
    return reads

# checks that the reads represent a continous region of the reference genome
def check_continous(reads):
    

def create_aligned_fasta(alignment_file):

    reads=extract_reads(alignment_file)
    for read in a.fetch():
        print(read.query_name)
        print(read.query_sequence)
        print(read.get_reference_sequence())


if __name__ == '__main__':
    # command line parser for the script. Current utilities are:
    # 1. create_aligned_fasta: create a fasta file from the aligned reads
    parser = argparse.ArgumentParser(description='Utilities for long read sequencing data analysis')
    parser.add_argument('utility', type=str, help="Utility to use,options are: " + ', '.join(utility_list) + '.')
    parser.add_argument('alignment_file', type=str, help='Alignment file')
    args = parser.parse_args()

    if args.utility == utility_list[0]:
        # log the utility being used to stderr
        sys.stderr.write('Using utility: ' + args.utility + '\n')

        create_aligned_fasta(args.alignment_file)