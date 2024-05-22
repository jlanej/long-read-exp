# utilities for long read sequencing data analysis
import os
import sys
import pysam
import argparse
import pandas as pd
from genomicranges import GenomicRanges

utility_list = ['create_aligned_fasta']


# extract all the reads from the alignment file to a list
def extract_reads(alignment_file):
    a = pysam.AlignmentFile(alignment_file, "rb")
    reads = []
    for read in a:
        reads.append(read)
    return reads


def create_pandas_df(reads):
    d = []
    for read in reads:
        d.append({'seqnames': read.reference_name, 'strand': read.is_reverse, 'cigar': read.cigarstring,
                  'read_name': read.query_name, 'starts': read.reference_start, 'ends': read.reference_end})

    return pd.DataFrame(d)


# checks that the reads represent a continous region of the reference genome
def convert_to_range(reads):
    gr = GenomicRanges.from_pandas(create_pandas_df(reads))
    #     ensure same seqnames
    if len(set(gr.get_seqnames())) > 1:
        sys.stderr.write('Error: Reads are not from the same chromosome\n')
        sys.exit(1)
    sys.stderr.write('Reads are from the same chromosome\n')
    print("prior to reduce", file=sys.stderr)
    print(gr, file=sys.stderr)
    gr_reduce = gr.reduce(ignore_strand=True)
    print("after reduce", file=sys.stderr)
    print(gr_reduce, file=sys.stderr)
    if len(gr_reduce) > 1:
        sys.stderr.write('Error: Reads span  discontinuous range\n')
        sys.exit(1)
    #     return the genomic range and the reduced genomic range
    return gr, gr_reduce


def create_aligned_fasta(alignment_file):
    reads = extract_reads(alignment_file)
    gr, gr_reduce = convert_to_range(reads)


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
