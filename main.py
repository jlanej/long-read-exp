# utilities for long read sequencing data analysis
import os
import sys
import pysam
import cigar
import argparse
import pandas as pd
from genomicranges import GenomicRanges
import matplotlib.pyplot as plt

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
    grr = GenomicRanges.from_pandas(create_pandas_df(reads))
    if len(set(grr.get_seqnames())) > 1:
        sys.stderr.write('Error: Reads are not from the same chromosome\n')
        sys.exit(1)
    sys.stderr.write('Reads are from the same chromosome\n')
    reduce = grr.reduce(ignore_strand=True)
    if len(reduce) > 1:
        sys.stderr.write('Error: Reads span  discontinuous range\n')
        sys.exit(1)
    return grr, reduce

# plot alignments as rectangles on a plot, one rectangle per read,
# with the read name on the y-axis and the rectangle spanning the start and end of the read

def plot_alignments(gr,gr_reduce):
    fig, ax = plt.subplots()
    for i in range(len(gr)):
        # ax.text(0, i, gr.mcols.get_column('read_name')[i], fontsize=12)
        ax.add_patch(plt.Rectangle((gr.get_start()[i], i - 0.4), gr.get_width()[i], 0.8, color='blue'))
    ax.set_xlim(gr_reduce.get_start()[0], gr_reduce.get_end()[0])
    ax.set_ylim(0, len(gr))
    plt.show()

def parse_haplotype(haplotype_file):
    reads = extract_reads(haplotype_file)
    gr, gr_reduce = convert_to_range(reads)
    return gr, gr_reduce


if __name__ == '__main__':
    # command line parser for the script. Current utilities are:
    # 1. create_aligned_fasta: create a fasta file from the aligned reads
    parser = argparse.ArgumentParser(description='Utilities for long read sequencing data analysis')
    parser.add_argument('utility', type=str, help="Utility to use,options are: " + ', '.join(utility_list) + '.')
    parser.add_argument('haplotype1', type=str, help='Alignment file for haplotype 1')
    parser.add_argument('haplotype2', type=str, help='Alignment file for haplotype 2')
    # output file
    parser.add_argument('output_file', type=str, help='Output file')
    args = parser.parse_args()

    if args.utility == utility_list[0]:
        # log the utility being used to stderr
        sys.stderr.write('Using utility: ' + args.utility + '\n')
        reads = extract_reads(args.alignment_file)
        gr, gr_reduce = convert_to_range(reads)
        plot_alignments(gr, gr_reduce)

