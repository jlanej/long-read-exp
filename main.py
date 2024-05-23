import sys

import cigar
import pysam
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


def get_non_clipped_cigar_length(cigar_string):
    #     Soft and Hard clipped bases are not included in the cigar string length
    cigar_s = cigar.Cigar(cigar_string)
    length = 0
    for c in cigar_s.items():
        if c[1] in ['M', 'I', 'X', '=']:
            length += c[0]
    return length


# for each read ID shared between the two haplotypes, choose the haplotype with longest alignment
def define_best_haplotype(gr1, gr2):
    gr1_ids = set(gr1.mcols.get_column('read_name'))
    # initialize a column to store the haplotype that the read belongs to
    gr1.mcols.set_column('belongs_to_this_hap', [False] * len(gr1), in_place=True)
    gr2_ids = set(gr2.mcols.get_column('read_name'))
    gr2.mcols.set_column('belongs_to_this_hap', [False] * len(gr2), in_place=True)
    common_ids = gr1_ids.intersection(gr2_ids)
    # map of read id to the best haplotype
    read_to_best_hap = {}

    if len(common_ids) == 0:
        sys.stderr.write('Error: No common read IDs between the two haplotypes\n')
        sys.exit(1)
    for read_id in common_ids:
        gr1_index = gr1.mcols.get_column('read_name').index(read_id)
        gr2_index = gr2.mcols.get_column('read_name').index(read_id)
        cigar1_len = get_non_clipped_cigar_length(gr1.mcols.get_column('cigar')[gr1_index])
        cigar2_len = get_non_clipped_cigar_length(gr2.mcols.get_column('cigar')[gr2_index])
        if cigar1_len > cigar2_len:
            gr1.mcols.get_column('belongs_to_this_hap')[gr1_index] = True
            read_to_best_hap[read_id] = 1
        elif cigar1_len < cigar2_len:
            gr2.mcols.get_column('belongs_to_this_hap')[gr2_index] = True
            read_to_best_hap[read_id] = 2
        else:
            #     set both to true if the cigar lengths are the same
            gr1.mcols.get_column('belongs_to_this_hap')[gr1_index] = True
            gr2.mcols.get_column('belongs_to_this_hap')[gr2_index] = True
            read_to_best_hap[read_id] = 0
    return gr1, gr2, read_to_best_hap


def plot_alignment(start, width, ax, color, y):
    ax.add_patch(plt.Rectangle((start, y - 0.4), width, 0.8, color=color))


# plot alignments as rectangles on a plot, one rectangle per read,
# with the read name on the y-axis and the rectangle spanning the start and end of the read

def plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original):
    fig, ax = plt.subplots(4)
    original_index = 0
    all_hap_index = 1
    best_hap_index = 2
    opposite_hap_index = 3

    for i in range(len(gr1)):
        if gr1.mcols.get_column('belongs_to_this_hap')[i]:
            plot_alignment(gr1.get_start()[i], gr1.get_width()[i], ax[best_hap_index], 'blue', i)
        else:
            plot_alignment(gr1.get_start()[i], gr1.get_width()[i], ax[opposite_hap_index], 'blue', i)
        plot_alignment(gr1.get_start()[i], gr1.get_width()[i], ax[all_hap_index], 'blue', i)
    for i in range(len(gr2)):
        if gr2.mcols.get_column('belongs_to_this_hap')[i]:
            plot_alignment(gr2.get_start()[i], gr2.get_width()[i], ax[best_hap_index], 'red', -1 * i)
        else:
            plot_alignment(gr2.get_start()[i], gr2.get_width()[i], ax[opposite_hap_index], 'red', -1 * i)
        plot_alignment(gr2.get_start()[i], gr2.get_width()[i], ax[all_hap_index], 'red', -1 * i)
    set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, best_hap_index)
    set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, all_hap_index)
    set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, opposite_hap_index)

    for i in range(len(gr_original)):
        plot_alignment(gr_original.get_start()[i], gr_original.get_width()[i], ax[original_index], 'green', i)
    ax[original_index].set_xlim(gr_reduce_original.get_start()[0], gr_reduce_original.get_end()[0])
    ax[original_index].set_ylim(0, len(gr_original))
    plt.show()


def set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, opposite_hap_index):
    ax[opposite_hap_index].set_xlim(min(gr_reduce1.get_start()[0], gr_reduce2.get_start()[0]),
                                    max(gr_reduce1.get_end()[0], gr_reduce2.get_end()[0]))
    ax[opposite_hap_index].set_ylim(-len(gr2), len(gr1))


def parse_haplotype(haplotype_file):
    reads = extract_reads(haplotype_file)
    gr, gr_reduce = convert_to_range(reads)
    return gr, gr_reduce, reads


if __name__ == '__main__':
    # command line parser for the script. Current utilities are:
    # 1. create_aligned_fasta: create a fasta file from the aligned reads
    parser = argparse.ArgumentParser(description='Utilities for long read sequencing data analysis')
    parser.add_argument('utility', type=str, help="Utility to use,options are: " + ', '.join(utility_list) + '.')
    parser.add_argument('haplotype1', type=str, help='Alignment file for haplotype 1')
    parser.add_argument('haplotype2', type=str, help='Alignment file for haplotype 2')
    parser.add_argument('original_alignment', type=str, help='Alignment file for original reads')

    args = parser.parse_args()

    if args.utility == utility_list[0]:
        # log the utility being used to stderr
        sys.stderr.write('Using utility: ' + args.utility + '\n')
        gr1, gr_reduce1, reads1 = parse_haplotype(args.haplotype1)
        gr2, gr_reduce2, reads2 = parse_haplotype(args.haplotype2)
        gr1, gr2, read_to_best_hap = define_best_haplotype(gr1, gr2)
        gr_original, gr_reduce_original, reads_original = parse_haplotype(args.original_alignment)
        plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original)
