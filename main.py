import os
import sys

import cigar
import pysam
import argparse
import pandas as pd
from genomicranges import GenomicRanges
import matplotlib.pyplot as plt

utility_list = ['create_aligned_fasta']


# https://github.com/moshi4/pyGenomeViz

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


def get_chr_start_stop(gr_row):
    return gr_row.get_seqnames()[0], gr_row.get_start()[0], gr_row.get_end()[0]


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


def get_non_clipped_cigar_span(cigar_string):
    #     Soft and Hard clipped bases are not included in the cigar string length
    cigar_s = cigar.Cigar(cigar_string)
    length = 0
    for c in cigar_s.items():
        if c[1] in ['M', 'I', 'X', '=']:
            length += c[0]
    return length


# find the indices of the first and last non-clipped bases in the read
def get_cigar_span_of_read(cigar_string):
    cigar_s = cigar.Cigar(cigar_string)
    start = 0
    end = 0
    for c in cigar_s.items():
        if c[1] in ['M', 'I', 'X', '=']:
            end += c[0]
        else:
            start += c[0]
    return start, end


def add_cigar_span(gr):
    cigar_span = []
    for i in range(len(gr)):
        cigar_span.append(get_non_clipped_cigar_span(gr.mcols.get_column('cigar')[i]))
    gr.mcols.set_column('cigar_span', cigar_span, in_place=True)
    return gr


def consolidate_cigar_spans_to_read_id(gr):
    #     get the cigar spans for each read
    gr = add_cigar_span(gr)
    #     get the cigar spans for each read
    read_name_to_cigar_span = {}
    read_name_ref_span = {}
    for i in range(len(gr)):
        read_name = gr.mcols.get_column('read_name')[i]
        if read_name not in read_name_to_cigar_span:
            read_name_to_cigar_span[read_name] = []
        if read_name not in read_name_ref_span:
            read_name_ref_span[read_name] = []
        read_name_to_cigar_span[read_name].append(gr.mcols.get_column('cigar_span')[i])
        read_name_ref_span[read_name].append(get_chr_start_stop(gr[i]))
    return read_name_to_cigar_span, read_name_ref_span


# transfer the span information from reads in gr1 to reads in gr2
def transfer_span_info(gr_source, gr_target):
    # #     add new columns to gr2 to store the span information
    # gr2.mcols.set_column('cigar_span_transfer', [0] * len(gr2), in_place=True)
    # gr2.mcols.set_column('start_transfer', [0] * len(gr2), in_place=True)
    # gr2.mcols.set_column('end_transfer', [0] * len(gr2),    in_place=True)

    # initiale the above columns to lists
    gr_target.mcols.set_column('cigar_span_transfer', [[]] * len(gr_target), in_place=True)
    gr_target.mcols.set_column('start_transfer', [[]] * len(gr_target), in_place=True)
    gr_target.mcols.set_column('end_transfer', [[]] * len(gr_target), in_place=True)

    read_name_to_cigar_span1, read_name_ref_span1 = consolidate_cigar_spans_to_read_id(gr_source)
    for read_name in read_name_to_cigar_span1:  # transfer the span information to gr2
        indices = [i for i, x in enumerate(gr_target.mcols.get_column('read_name')) if x == read_name]
        # print(len(indices))
        # print(indices)
        # indices = gr2.mcols.get_column('read_name').index(read_name)
        # print(indices)
        # print(len(indices))
        for i in indices:
            # print(gr2.mcols.get_column('start_transfer'))
            # print(len(gr2.mcols.get_column('start_transfer')))
            # print(gr2.mcols.get_column('start_transfer')[i])
            # print( read_name_ref_span1[read_name])
            # # print the type of read_name_ref_span1[read_name]
            # print(type(read_name_ref_span1[read_name]))
            # print(len(read_name_ref_span1[read_name]))

            # append all the values in read_name_ref_span1[read_name]
            for ref_span in read_name_ref_span1[read_name]:
                gr_target.mcols.get_column('cigar_span_transfer')[i].append(read_name_to_cigar_span1[read_name])
                gr_target.mcols.get_column('start_transfer')[i].append(ref_span[1])
                gr_target.mcols.get_column('end_transfer')[i].append(ref_span[2])

            # gr2.mcols.get_column('cigar_span_transfer')[i] = read_name_to_cigar_span1[read_name]
            # gr2.mcols.get_column('start_transfer')[i] = read_name_ref_span1[read_name][0][1]
            # gr2.mcols.get_column('end_transfer')[i] = read_name_ref_span1[read_name][0][2]
    return gr_target


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
        cigar1_len = get_non_clipped_cigar_span(gr1.mcols.get_column('cigar')[gr1_index])
        cigar2_len = get_non_clipped_cigar_span(gr2.mcols.get_column('cigar')[gr2_index])
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

def plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original, output_root):
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
    # plt.savefig('haplotype_alignment.png')
    print("saving to", output_root + '_haplotype_alignment.png')
    plt.savefig(output_root + '_haplotype_alignment.png', dpi=300, bbox_inches='tight')


def set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, opposite_hap_index):
    ax[opposite_hap_index].set_xlim(min(gr_reduce1.get_start()[0], gr_reduce2.get_start()[0]),
                                    max(gr_reduce1.get_end()[0], gr_reduce2.get_end()[0]))
    ax[opposite_hap_index].set_ylim(-len(gr2), len(gr1))


def parse_haplotype(haplotype_file):
    reads = extract_reads(haplotype_file)
    gr, gr_reduce = convert_to_range(reads)
    return gr, gr_reduce, reads


def write_reads_to_best_haplotype(reads, read_to_best_hap, output_root, header, write_both_haplotypes=False):
    h1_out_bam = output_root + '_hap1.bam'
    h2_out_bam = output_root + '_hap2.bam'
    h1_out = pysam.AlignmentFile(h1_out_bam, "wb", header=header)
    h2_out = pysam.AlignmentFile(h2_out_bam, "wb", header=header)
    for read in reads:
        if read.query_name in read_to_best_hap:
            if read_to_best_hap[read.query_name] == 1:
                h1_out.write(read)
            elif read_to_best_hap[read.query_name] == 2:
                h2_out.write(read)
            elif read_to_best_hap[read.query_name] == 0 and write_both_haplotypes:
                h1_out.write(read)
                h2_out.write(read)
    # close the output bam files
    h1_out.close()
    h2_out.close()
    # index the output bam files
    pysam.index(h1_out_bam)
    pysam.index(h2_out_bam)


def get_header_from_bam(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    return bam.header


if __name__ == '__main__':
    # command line parser for the script. Current utilities are:
    # 1. create_aligned_fasta: create a fasta file from the aligned reads
    parser = argparse.ArgumentParser(description='Utilities for long read sequencing data analysis')
    parser.add_argument('utility', type=str, help="Utility to use,options are: " + ', '.join(utility_list) + '.')
    parser.add_argument('haplotype1', type=str, help='Alignment file for haplotype 1')
    parser.add_argument('haplotype2', type=str, help='Alignment file for haplotype 2')
    parser.add_argument('original_alignment', type=str, help='Alignment file for original reads')
    # output directory argument
    parser.add_argument('output_directory', type=str, help='Directory to write output files to')
    args = parser.parse_args()

    if args.utility == utility_list[0]:
        # log the utility being used to stderr
        sys.stderr.write('Using utility: ' + args.utility + '\n')
        gr1, gr_reduce1, reads1 = parse_haplotype(args.haplotype1)
        gr2, gr_reduce2, reads2 = parse_haplotype(args.haplotype2)
        gr1, gr2, read_to_best_hap = define_best_haplotype(gr1, gr2)
        gr_original, gr_reduce_original, reads_original = parse_haplotype(args.original_alignment)

        root_original = args.output_directory + os.path.basename(args.original_alignment)
        transfer_span_info(gr_original, gr1)
        transfer_span_info(gr_original, gr2)
        plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original, root_original)
        #         output directory should be the same as the original alignment file
        #         output_directory = os.path.dirname(args.original_alignment)

        # root Should be the basename of the alignment file
        print("writing to root_original", root_original)
        write_reads_to_best_haplotype(reads_original, read_to_best_hap, root_original,
                                      get_header_from_bam(args.original_alignment))

        root_haplotype1 = args.output_directory + os.path.basename(args.haplotype1)
        print("writing to root_haplotype1", root_haplotype1)
        write_reads_to_best_haplotype(reads1, read_to_best_hap, root_haplotype1,
                                      get_header_from_bam(args.haplotype1))
        root_haplotype2 = args.output_directory + os.path.basename(args.haplotype2)
        print("writing to root_haplotype2", root_haplotype2)
        write_reads_to_best_haplotype(reads2, read_to_best_hap, root_haplotype2,
                                      get_header_from_bam(args.haplotype2))
        # plot_haplotypes_plotly(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original, root_original)
