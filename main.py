import os
import sys

import cigar
import pysam
import argparse
import matplotlib.pyplot as plt

import utils.hap_cluster as hap_cluster
import utils.lr_utils as lr_utils



utility_list = ['create_aligned_fasta']


# https://github.com/moshi4/pyGenomeViz

# extract all the reads from the alignment file to a list

def get_cigar_length_for_ops(cigar_string):
    subtract_ops = ['D', 'I', 'X', 'S', 'H']
    ops = ['M', '=']
    #     Soft and Hard clipped bases are not included in the cigar string length
    cigar_s = cigar.Cigar(cigar_string)
    length = 0
    for c in cigar_s.items():
        if c[1] in ops:
            length += c[0]
        elif c[1] in subtract_ops:
            length -= c[0]
        else:
            sys.stderr.write('Error: Unknown cigar operation\n' + c[1] + '\n')
            sys.exit(1)
    return length



def reverse_span(span, length):
    return abs(span[0] - length), abs(span[1] - length)


# for each read ID shared between the two haplotypes, choose the haplotype with longest alignment
def define_best_haplotype(gr1, gr2):
    gr1_ids = set(gr1.mcols.get_column('read_name'))
    # initialize a column to store the haplotype that the read belongs to
    gr1.mcols.set_column('belongs_to_this_hap', [False] * len(gr1), in_place=True)
    gr2_ids = set(gr2.mcols.get_column('read_name'))
    gr2.mcols.set_column('belongs_to_this_hap', [False] * len(gr2), in_place=True)
    common_ids = gr1_ids.intersection(gr2_ids)
    # warn about reads that are not in both haplotypes
    lr_utils.warn_diff_read_ids(gr1_ids, gr2_ids)

    # map of read id to the best haplotype
    read_to_best_hap = {}

    if len(common_ids) == 0:
        sys.stderr.write('Error: No common read IDs between the two haplotypes\n')
        sys.exit(1)
    for read_id in common_ids:
        gr1_index = gr1.mcols.get_column('read_name').index(read_id)
        gr2_index = gr2.mcols.get_column('read_name').index(read_id)
        cigar1_len = get_cigar_length_for_ops(gr1.mcols.get_column('cigar')[gr1_index])
        cigar2_len = get_cigar_length_for_ops(gr2.mcols.get_column('cigar')[gr2_index])
        if cigar1_len > cigar2_len:
            gr1.mcols.get_column('belongs_to_this_hap')[gr1_index] = True
            read_to_best_hap[read_id] = 1
        elif cigar1_len < cigar2_len:
            gr2.mcols.get_column('belongs_to_this_hap')[gr2_index] = True
            read_to_best_hap[read_id] = 2
        else:
            # set both to true if the cigar lengths are the same
            gr1.mcols.get_column('belongs_to_this_hap')[gr1_index] = True
            gr2.mcols.get_column('belongs_to_this_hap')[gr2_index] = True
            read_to_best_hap[read_id] = 0
    return gr1, gr2, read_to_best_hap



def plot_alignment(start, width, ax, color, y):
    ax.add_patch(plt.Rectangle((start, y - 0.4), width, 0.8, color=color))


def get_n_colors(n):
    return plt.colormaps['viridis'].resampled(n).colors


# plots each portion of the reads's alignment span as a different colred rectangle
def plot_alignment_span(gr_row, ax, y, read_name_to_cigar_span, read_name_ref_span):
    read_name = gr_row.mcols.get_column('read_name')[0]
    spans = read_name_to_cigar_span[read_name]

    spans = sorted(spans, key=lambda x: x[0])
    colors = get_n_colors(len(spans))

    if len(spans) > 0:
        for i in range(len(spans)):
            start, end = spans[i]
            start = start + gr_row.get_start()[0]
            end = end + gr_row.get_start()[0]
            if gr_row.get_strand()[0] == 1:
                ax.add_patch(plt.Rectangle((start, y - 0.4), end - start, 0.4, color=colors[i], edgecolor='black'))
            else:
                ax.add_patch(plt.Rectangle((start, y - 0.4), end - start, 0.4, color=colors[i], edgecolor='black'))


def plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original, output_root):
    fig, ax = plt.subplots(5)
    test_index = 4
    test_index_opp = 5

    original_index = 0
    all_hap_index = 1
    best_hap_index = 3
    opposite_hap_index = 2

    read_name_to_cigar_span, read_name_ref_span = lr_utils.consolidate_cigar_spans_to_read_id(gr_original)

    h1_plot_index = 0
    h2_plot_index = 0
    for i in range(len(gr1)):
        if gr1.mcols.get_column('belongs_to_this_hap')[i]:
            plot_alignment(gr1.get_start()[i], gr1.get_width()[i], ax[best_hap_index], 'blue', h1_plot_index)
            plot_alignment_span(gr1[i], ax[test_index], h1_plot_index, read_name_to_cigar_span, read_name_ref_span)
            h1_plot_index += 1
        else:
            plot_alignment(gr1.get_start()[i], gr1.get_width()[i], ax[opposite_hap_index], 'blue', i)
            # plot_alignment_span(gr1[i], ax[test_index_opp], i, read_name_to_cigar_span, read_name_ref_span)
        plot_alignment(gr1.get_start()[i], gr1.get_width()[i], ax[all_hap_index], 'blue', i)
    for i in range(len(gr2)):
        if gr2.mcols.get_column('belongs_to_this_hap')[i]:
            plot_alignment(gr2.get_start()[i], gr2.get_width()[i], ax[best_hap_index], 'red', -1 * h2_plot_index)
            plot_alignment_span(gr2[i], ax[test_index], -1 * h2_plot_index, read_name_to_cigar_span, read_name_ref_span)
            h2_plot_index += 1
        else:
            plot_alignment(gr2.get_start()[i], gr2.get_width()[i], ax[opposite_hap_index], 'red', -1 * i)
            # plot_alignment_span(gr2[i], ax[test_index_opp], -1 * i, read_name_to_cigar_span, read_name_ref_span)
        plot_alignment(gr2.get_start()[i], gr2.get_width()[i], ax[all_hap_index], 'red', -1 * i)
    set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, best_hap_index)
    set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, all_hap_index)
    set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, opposite_hap_index)
    set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, test_index)
    # set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, test_index_opp)
    ax[test_index].set_ylim(-1 * (5 + h2_plot_index), h1_plot_index + 5)
    ax[best_hap_index].set_ylim(-1 * (h2_plot_index + 5), h1_plot_index + 5)

    for i in range(len(gr_original)):
        plot_alignment(gr_original.get_start()[i], gr_original.get_width()[i], ax[original_index], 'green', i)
    ax[original_index].set_xlim(gr_reduce_original.get_start()[0], gr_reduce_original.get_end()[0])
    ax[original_index].set_ylim(0, len(gr_original))
    # plt.savefig('haplotype_alignment.png')
    print("saving to", output_root + '_haplotype_alignment.png')
    plt.gcf().set_size_inches(10, 15)
    plt.savefig(output_root + '_haplotype_alignment.png', dpi=300, bbox_inches='tight')


def set_axis(ax, gr1, gr2, gr_reduce1, gr_reduce2, opposite_hap_index):
    ax[opposite_hap_index].set_xlim(min(gr_reduce1.get_start()[0], gr_reduce2.get_start()[0]),
                                    max(gr_reduce1.get_end()[0], gr_reduce2.get_end()[0]))
    ax[opposite_hap_index].set_ylim(-len(gr2), len(gr1))


def parse_haplotype(haplotype_file):
    reads = lr_utils.extract_reads(haplotype_file)
    gr, gr_reduce = lr_utils.convert_to_range(reads)
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
    hap_cluster.hello()
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
        plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original, root_original)
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
