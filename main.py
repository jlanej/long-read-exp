import os
import sys

import cigar
import pysam
import argparse
import matplotlib.pyplot as plt
import utils.hap_cluster as hap_cluster
import utils.lr_utils as lr_utils

utility_list = ['create_aligned_fasta']


# for each read ID shared between the two haplotypes, choose the haplotype with longest alignment
def define_best_haplotype(gr1, gr2):
    gr1_ids = set(gr1.mcols.get_column('read_name'))
    gr2_ids = set(gr2.mcols.get_column('read_name'))
    common_ids = gr1_ids.intersection(gr2_ids)
    lr_utils.warn_diff_read_ids(gr1_ids, gr2_ids, exit_on_diff=True)
    # map of read id to the best haplotype
    read_to_best_hap = {}

    if len(common_ids) == 0:
        sys.stderr.write('Error: No common read IDs between the two haplotypes\n')
        sys.exit(1)
    for read_id in common_ids:
        gr1_index = gr1.mcols.get_column('read_name').index(read_id)
        gr2_index = gr2.mcols.get_column('read_name').index(read_id)
        cigar1_len = lr_utils.get_cigar_length_metric(gr1.mcols.get_column('cigar')[gr1_index])
        cigar2_len = lr_utils.get_cigar_length_metric(gr2.mcols.get_column('cigar')[gr2_index])
        if cigar1_len > cigar2_len:
            read_to_best_hap[read_id] = 1
        elif cigar1_len < cigar2_len:
            read_to_best_hap[read_id] = 2
        else:
            read_to_best_hap[read_id] = 0
    return gr1, gr2, read_to_best_hap


def plot_alignment(start, width, ax, color, y):
    ax.add_patch(plt.Rectangle((start, y - 0.25), width, 0.7, color=color))


def get_n_viridis_colors(n, pal):
    return plt.colormaps[pal].resampled(n).colors


# plots each portion of the reads's alignment span as a different colred rectangle
def plot_alignment_span(gr_row, ax, y, read_name_to_cigar_span, pal):
    read_name = gr_row.mcols.get_column('read_name')[0]
    spans = read_name_to_cigar_span[read_name]
    spans = sorted(spans, key=lambda x: x[0])
    colors = get_n_viridis_colors(len(spans) + 1, pal)

    for i in range(len(spans)):
        start, end = spans[i]
        start = start + gr_row.get_start()[0]
        end = end + gr_row.get_start()[0]
        # if the length of the span is more than 1, the color is i+1
        color = colors[i]
        if len(spans) > 1:
            color = colors[i + 1]
        plot_alignment(start, end - start, ax, color, y)


def plot_reference_span(gr_row, ax, y, read_name_to_cigar_span, read_name_ref_span, pal):
    read_name = gr_row.mcols.get_column('read_name')[0]
    spans = read_name_to_cigar_span[read_name]
    ref_spans = read_name_ref_span[read_name]
    colors = get_n_viridis_colors(len(spans) + 1, pal)

    ref_sort = [x for _, x in sorted(zip(spans, ref_spans), key=lambda pair: pair[0][0])]
    spans = sorted(spans, key=lambda x: x[0])
    cigar_ops = cigar.Cigar(gr_row.mcols.get_column('cigar')[0])
    span_current = lr_utils.get_cigar_span_of_read(cigar_ops)
    span_index = spans.index(span_current)
    chrom, start, end = ref_sort[span_index]
    color = colors[span_index]
    if len(spans) > 1:
        color = colors[span_index + 1]
    plot_alignment(start, end - start, ax, color, y)


def plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original, output_root,
                    read_name_to_cigar_metrics):
    fig, ax = plt.subplots(3)

    original_index = 0
    assembly_index_h1 = 1
    assembly_index_h2 = 2

    read_name_to_cigar_span, read_name_ref_span = lr_utils.consolidate_cigar_spans_to_read_id(gr_original)
    h1_pal = 'magma'
    h2_pal = 'viridis'
    index_start = 5
    h1_plot_index = index_start
    h2_plot_index = index_start

    for i in range(len(gr1)):
        if lr_utils.get_best_hap_for_read(gr1.mcols.get_column('read_name')[i], read_name_to_cigar_metrics) == 1:
            plot_alignment_span(gr1[i], ax[assembly_index_h1], h1_plot_index, read_name_to_cigar_span, h1_pal)
            h1_plot_index += 1
    for i in range(len(gr2)):
        if lr_utils.get_best_hap_for_read(gr2.mcols.get_column('read_name')[i], read_name_to_cigar_metrics) == 2:
            plot_alignment_span(gr2[i], ax[assembly_index_h2], -1 * h2_plot_index, read_name_to_cigar_span, h2_pal)
            h2_plot_index += 1

    set_axis(ax[assembly_index_h1], 0, h1_plot_index + index_start,gr_reduce1.get_start()[0], gr_reduce1.get_end()[0], 'Haplotype 1 Alignment')
    label_hap(ax[assembly_index_h1], h1_plot_index/2, 'H1')
    set_axis(ax[assembly_index_h2], -1 * (index_start + h2_plot_index), 0,
             gr_reduce2.get_start()[0], gr_reduce2.get_end()[0], 'Haplotype 2 Alignment')
    label_hap(ax[assembly_index_h2], -1 * h2_plot_index/2, 'H2')

    h1_plot_index = index_start
    h2_plot_index = index_start
    for i in range(len(gr_original)):
        if lr_utils.get_best_hap_for_read(gr_original.mcols.get_column('read_name')[i],
                                          read_name_to_cigar_metrics) == 1:
            plot_reference_span(gr_original[i], ax[original_index], h1_plot_index, read_name_to_cigar_span,
                                read_name_ref_span, h1_pal)
            h1_plot_index += 1
        elif lr_utils.get_best_hap_for_read(gr_original.mcols.get_column('read_name')[i],
                                            read_name_to_cigar_metrics) == 2:
            plot_reference_span(gr_original[i], ax[original_index], -1 * h2_plot_index, read_name_to_cigar_span,
                                read_name_ref_span, h2_pal)
            h2_plot_index += 1

    set_axis(ax[original_index], -1 * (index_start + h2_plot_index), h1_plot_index + index_start,
             min(gr_reduce_original.get_start()[0], gr_reduce_original.get_start()[0]),
             max(gr_reduce_original.get_end()[0], gr_reduce_original.get_end()[0]),
             'Original Reference Alignment')
    label_haps(ax[original_index], h2_plot_index, h1_plot_index)

    print("saving to", output_root + '_haplotype_alignment.png')
    plt.gcf().set_size_inches(15, 15)
    plt.savefig(output_root + '_haplotype_alignment.png', dpi=300, bbox_inches='tight')


def set_axis(ax, ymin, ymax, xmin, xmax, title):
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axhline(y=0, color='grey', linestyle='dotted', linewidth=4)
    ax.set_title(title)
def label_hap(ax,index, hap):
    ax.set_yticks([index])
    ax.set_yticklabels([hap])

def label_haps(ax, h1_index, h2_index):
    ax.set_yticks([-1 * h2_index / 2, h1_index / 2])
    ax.set_yticklabels(['H2', 'H1'])


def parse_haplotype(haplotype_file):
    reads = lr_utils.extract_reads(haplotype_file)
    gr, gr_reduce = lr_utils.convert_to_range(reads)
    return gr, gr_reduce, reads


def write_reads_to_best_haplotype(reads, read_name_to_cigar_metrics, output_root, header, write_both_haplotypes=False):
    h1_out_bam = output_root + '_hap1.bam'
    h2_out_bam = output_root + '_hap2.bam'
    h_combined_out_bam = output_root + '_hap_combined.bam'
    # prep the header for two new read tags : best haplotype and cigar metrics
    h1_out = pysam.AlignmentFile(h1_out_bam, "wb", header=header)
    h2_out = pysam.AlignmentFile(h2_out_bam, "wb", header=header)
    h_combined_out = pysam.AlignmentFile(h_combined_out_bam, "wb", header=header)

    for read in reads:
        read.set_tag('BH', lr_utils.get_best_hap_for_read(read.query_name, read_name_to_cigar_metrics))
        read.set_tag('C1', read_name_to_cigar_metrics.loc[read.query_name, "hap1_prop"], value_type='f')
        read.set_tag('C2', read_name_to_cigar_metrics.loc[read.query_name, "hap2_prop"], value_type='f')
        read.set_tag('M1', read_name_to_cigar_metrics.loc[read.query_name, "hap1"], value_type='f')
        read.set_tag('M2', read_name_to_cigar_metrics.loc[read.query_name, "hap2"], value_type='f')

        if read.query_name in read_name_to_cigar_metrics.index:
            h_combined_out.write(read)
            if lr_utils.get_best_hap_for_read(read.query_name, read_name_to_cigar_metrics) == 1:
                h1_out.write(read)
            elif lr_utils.get_best_hap_for_read(read.query_name, read_name_to_cigar_metrics) == 2:
                h2_out.write(read)
            elif lr_utils.get_best_hap_for_read(read.query_name,
                                                read_name_to_cigar_metrics) == 0 and write_both_haplotypes:
                h1_out.write(read)
                h2_out.write(read)
    # close the output bam files
    h1_out.close()
    h2_out.close()
    h_combined_out.close()
    # index the output bam files
    pysam.index(h1_out_bam)
    pysam.index(h2_out_bam)
    pysam.index(h_combined_out_bam)


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
    # paf file argument
    parser.add_argument('--paf', type=str, help='PAF file to read v read alignment information')
    args = parser.parse_args()

    if args.utility == utility_list[0]:
        # log the utility being used to stderr
        sys.stderr.write('Using utility: ' + args.utility + '\n')
        # create the output directory if it does not exist
        if not os.path.exists(args.output_directory):
            os.makedirs(args.output_directory)

        gr1, gr_reduce1, reads1 = parse_haplotype(args.haplotype1)
        gr2, gr_reduce2, reads2 = parse_haplotype(args.haplotype2)

        # indels smaller than this size will be ignored for the purposes of haplotype assignment
        min_indel_size = 2
        read_name_to_cigar_metrics = lr_utils.prep_hap_metrics(gr1, gr2, min_indel_size)

        gr_original, gr_reduce_original, reads_original = parse_haplotype(args.original_alignment)
        root_original = args.output_directory + os.path.basename(args.original_alignment)
        plot_haplotypes(gr1, gr2, gr_reduce1, gr_reduce2, gr_original, gr_reduce_original, root_original,
                        read_name_to_cigar_metrics)
        print("writing to root_original", root_original)
        write_reads_to_best_haplotype(reads_original, read_name_to_cigar_metrics, root_original,
                                      get_header_from_bam(args.original_alignment))
        #
        root_haplotype1 = args.output_directory + os.path.basename(args.haplotype1)
        print("writing to root_haplotype1", root_haplotype1)
        write_reads_to_best_haplotype(reads1, read_name_to_cigar_metrics, root_haplotype1,
                                      get_header_from_bam(args.haplotype1))
        root_haplotype2 = args.output_directory + os.path.basename(args.haplotype2)
        print("writing to root_haplotype2", root_haplotype2)
        write_reads_to_best_haplotype(reads2, read_name_to_cigar_metrics, root_haplotype2,
                                      get_header_from_bam(args.haplotype2))
        lr_utils.cluster_haplotypes(gr1, gr2, read_name_to_cigar_metrics, root_original + '_haplotype_cluster.png')
