# utilities for long read sequencing data analysis
import os
import sys
import pysam
import cigar
import argparse
import pandas as pd
from genomicranges import GenomicRanges
import dash_bio as dashbio
from dash import Dash, html, Input, Output, callback
import urllib.request as urlreq

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
        seq_name = read.reference_name
        # if the sequence name is a non-standard chromosome, add chr0 to the sequence name There appears to be a bug
        # in the reduce function that causes it to fail when the sequence name is not a standard chromosome
        if not seq_name.startswith('chr'):
            # replace "_" with "." in the sequence name
            seq_name = seq_name.replace("_", ".")
            seq_name = seq_name + '_0'
        d.append({'seqnames': seq_name, 'strand': read.is_reverse, 'cigar': read.cigarstring,
                  'read_name': read.query_name, 'starts': read.reference_start, 'ends': read.reference_end,
                  'sequence': read.query_sequence})

    return pd.DataFrame(d)


# checks that the reads represent a continous region of the reference genome
def convert_to_range(reads):
    gr = GenomicRanges.from_pandas(create_pandas_df(reads))
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
    return gr, gr_reduce


# makes all read sequences the same length by adding Ns to the beginning or end of the sequence

def harmonize_sequences(gr, gr_reduce):
    all_start = gr_reduce.get_start()[0]
    all_end = gr_reduce.get_end()[0]
    gr.mcols.set_column('original_sequence', gr.mcols.get_column('sequence'), in_place=True)

    for i in range(len(gr)):
        cigars = cigar.Cigar(gr.mcols.get_column('cigar')[i])
        start = gr.get_start()[i]
        current_sequence = gr.mcols.get_column('sequence')[i]
        num_Ns = start - all_start
        if start > all_start:
            current_sequence = '.' * num_Ns + current_sequence

        indexToInsert = num_Ns
        for c in cigars.items():
            #     if the operation consumes the reference genome but not the read, add "." to the read sequence
            #     to match the reference genome
            if c[1] in ['D']:
                current_sequence = current_sequence[:indexToInsert] + '.' * c[0] + current_sequence[indexToInsert:]
                indexToInsert += c[0]
        effective_end = len(current_sequence)

        if effective_end < all_end:
            current_sequence = current_sequence + '.' * (all_end - effective_end)
        if len(current_sequence) < all_end - all_start:
            sys.stderr.write('Error: Sequence length too short\n')
            # print the sequence length and the expected length
            print(effective_end)
            print(len(current_sequence), all_end - all_start, file=sys.stderr)
            sys.exit(1)
        # trim the sequence to 10000 bases on either side of the middle of the full range
        middle = (all_end - all_start) // 2
        # current_sequence = current_sequence[middle - 10000:middle + 10000]
        gr.mcols.get_column('sequence')[i] = current_sequence

    return gr


def create_harmonized_sequence(alignment_file):
    reads = extract_reads(alignment_file)
    gr, gr_reduce = convert_to_range(reads)
    gr = harmonize_sequences(gr, gr_reduce)
    return gr


def write_harmonized_fasta(gr, output_file):
    out = open(output_file, 'w')
    for i in range(len(gr)):
        out.write('>' + gr.get_seqnames()[0] + ':' + str(gr.get_start()[0]) + '-' + str(gr.get_end()[0]) + ':' +
                  str(gr.get_strand()[i]) + ':' + gr.mcols.get_column('read_name')[i] + '\n')
        out.write(gr.mcols.get_column('sequence')[i] + '\n')
    out.close()


if __name__ == '__main__':
    # command line parser for the script. Current utilities are:
    # 1. create_aligned_fasta: create a fasta file from the aligned reads
    parser = argparse.ArgumentParser(description='Utilities for long read sequencing data analysis')
    parser.add_argument('utility', type=str, help="Utility to use,options are: " + ', '.join(utility_list) + '.')
    parser.add_argument('alignment_file', type=str, help='Alignment file')
    # output file
    parser.add_argument('output_file', type=str, help='Output file')
    args = parser.parse_args()

    if args.utility == utility_list[0]:
        # log the utility being used to stderr
        sys.stderr.write('Using utility: ' + args.utility + '\n')
        gr = create_harmonized_sequence(args.alignment_file)
        write_harmonized_fasta(gr, args.output_file)

