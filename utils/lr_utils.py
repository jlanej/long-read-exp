import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pysam
from genomicranges import GenomicRanges
import cigar
from collections import Counter


def warn_diff_read_ids(gr1_ids, gr2_ids, exit_on_diff=True):
    if len(gr1_ids) != len(gr2_ids):
        sys.stderr.write('Warning: Reads are not the same between the two haplotypes\n')
        #     list the read IDs that are not in both haplotypes
        for read_id in gr1_ids.difference(gr2_ids):
            sys.stderr.write('Read ID: ' + read_id + ' not in haplotype 2\n')
            if exit_on_diff:
                sys.exit(1)
        for read_id in gr2_ids.difference(gr1_ids):
            sys.stderr.write('Read ID: ' + read_id + ' not in haplotype 1\n')
            if exit_on_diff:
                sys.exit(1)


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


# checks that the reads represent a continuous region of the reference genome
def convert_to_range(reads, keep_most_common_chr=True, exit_on_error=True):
    grr = GenomicRanges.from_pandas(create_pandas_df(reads))
    if len(set(grr.get_seqnames())) > 1:
        sys.stderr.write('Error: Reads are not from the same chromosome\n')
        if keep_most_common_chr:
            chr_counts = Counter(grr.get_seqnames())
            chr_counts = dict(chr_counts)
            most_common_chr = max(chr_counts, key=chr_counts.get)
            indices_to_keep = [i for i, x in enumerate(grr.get_seqnames()) if x == most_common_chr]
            grr = grr.get_subset(indices_to_keep)
            sys.stderr.write('Keeping reads from the most common chromosome: ' + most_common_chr + '\n')
        elif exit_on_error:
            sys.exit(1)
    reduce = grr.reduce(ignore_strand=True)
    if len(reduce) > 1:
        sys.stderr.write('Error: Reads span  discontinuous range\n')
        if exit_on_error:
            sys.exit(1)
    return grr, reduce


def add_cigar_span(gr):
    cigar_span = []
    for i in range(len(gr)):
        cigar_ops = cigar.Cigar(gr.mcols.get_column('cigar')[i])
        span = get_cigar_span_of_read(cigar_ops)
        # if gr.get_strand()[i] == 1:
        #    span = reverse_span(span, gr.get_width()[i])
        cigar_span.append(span)
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


# find the indices of the first and last non-clipped bases in the read
def get_cigar_span_of_read(cigar_s):
    start = 0
    end = 0
    finding_start = True
    for c in cigar_s.items():
        if c[1] in ['M', 'I', 'X', '=']:
            finding_start = False
            end += c[0]
        elif c[1] in ['S', 'H'] and finding_start:
            start += c[0]
    return start, end + start


def reverse_span(span, length):
    return abs(span[0] - length), abs(span[1] - length)


# IGV hides indels with length less than 2, and since these are common and noisy, we will ignore them for now as well
def get_cigar_length_metric(cigar_string, min_indel_size):
    subtract_ops = ['D', 'I', 'X', 'S', 'H']
    ops = ['M', '=']
    cigar_s = cigar.Cigar(cigar_string)
    length = 0
    for c in cigar_s.items():
        if c[1] in ops:
            length += c[0]
        elif c[1] in subtract_ops:
            if (c[1] != 'I' and c[1] != 'D') or c[0] >= min_indel_size:
                length -= c[0]
        else:
            sys.stderr.write('Error: Unknown cigar operation\n' + c[1] + '\n')
            sys.exit(1)
    return length


def get_cigar_match_metric(cigar_string):
    ops = ['M', '=']
    cigar_s = cigar.Cigar(cigar_string)
    length = 0
    for c in cigar_s.items():
        if c[1] in ops:
            length += c[0]
    return length


def get_spans_per_read(gr, min_indel_size):
    read_name_to_length_metric = {}
    read_name_to_match_metric = {}
    for i in range(len(gr)):
        read_id = gr.mcols.get_column('read_name')[i]
        if read_id not in read_name_to_length_metric:
            read_name_to_length_metric[read_id] = 0
        if read_id not in read_name_to_match_metric:
            read_name_to_match_metric[read_id] = 0
        read_name_to_length_metric[read_id] += get_cigar_length_metric(gr.mcols.get_column('cigar')[i], min_indel_size)
        read_name_to_match_metric[read_id] += get_cigar_match_metric(gr.mcols.get_column('cigar')[i])
    return read_name_to_length_metric, read_name_to_match_metric


def get_read_length(gr):
    read_name_to_seq_len = {}
    for i in range(len(gr)):
        read_id = gr.mcols.get_column('read_name')[i]
        read_name_to_seq_len[read_id] = get_full_read_length_from_cigar(gr.mcols.get_column('cigar')[i])
    return read_name_to_seq_len


def get_full_read_length_from_cigar(cigar_string):
    cigar_s = cigar.Cigar(cigar_string)
    length = 0
    for c in cigar_s.items():
        if c[1] in ['M', 'I', 'X', '=', 'S', 'H']:
            length += c[0]
    return length


# Col	Type	Description
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start coordinate (0-based)
# 4	int	Query end coordinate (0-based)
# 5	char	‘+’ if query/target on the same strand; ‘-’ if opposite
# 6	string	Target sequence name
# 7	int	Target sequence length
# 8	int	Target start coordinate on the original strand
# 9	int	Target end coordinate on the original strand
# 10	int	Number of matching bases in the mapping
# 11	int	Number bases, including gaps, in the mapping
# 12	int	Mapping quality (0-255 with 255 for missing)

def load_paf_file(paf_file):
    paf = pd.read_csv(paf_file, sep='\t', header=None)
    paf.columns = ['query_name', 'query_length', 'query_start', 'query_end', 'strand', 'target_name',
                   'target_length', 'target_start', 'target_end', 'num_matches', 'num_bases', 'mapping_quality']
    return paf


def prep_hap(gr, hap_num, min_indel_size):
    read_name_to_cigar_span, read_name_to_match_metric = get_spans_per_read(gr, min_indel_size)
    read_name_to_cigar_span = pd.DataFrame.from_dict(read_name_to_cigar_span, orient='index')
    read_name_to_cigar_span.columns = ['hap' + str(hap_num)]
    read_name_to_match_metric = pd.DataFrame.from_dict(read_name_to_match_metric, orient='index')
    read_name_to_match_metric.columns = ['hap' + str(hap_num) + '_match_metric']
    read_name_to_cigar_span = read_name_to_cigar_span.merge(read_name_to_match_metric, left_index=True,
                                                            right_index=True, how='outer')
    return read_name_to_cigar_span


def prep_hap_metrics(gr1, gr2, min_indel_size):
    read_name_to_cigar_span1 = prep_hap(gr1, 1, min_indel_size)
    read_name_to_cigar_span2 = prep_hap(gr2, 2, min_indel_size)
    warn_diff_read_ids(set(read_name_to_cigar_span1.keys()), set(read_name_to_cigar_span2.keys()))

    # merge the two data frames on the read name
    read_name_to_cigar_metrics = read_name_to_cigar_span1.merge(read_name_to_cigar_span2, left_index=True,
                                                                right_index=True, how='outer')

    read_name_seq_len1 = pd.DataFrame.from_dict(get_read_length(gr1), orient='index')
    read_name_seq_len1.columns = ['read_length']
    # merge the two data frames on the read name
    read_name_to_cigar_metrics = read_name_to_cigar_metrics.merge(read_name_seq_len1, left_index=True, right_index=True,
                                                                  how='outer')

    props_for_proportion = ['hap1', 'hap2', 'hap1_match_metric', 'hap2_match_metric']
    proportions = ['hap1_prop', 'hap2_prop', 'hap1_match_prop', 'hap2_match_prop']

    for i in range(len(props_for_proportion)):
        read_name_to_cigar_metrics[proportions[i]] = read_name_to_cigar_metrics[props_for_proportion[i]] / \
                                                     read_name_to_cigar_metrics['read_length']
        read_name_to_cigar_metrics[proportions[i]] = read_name_to_cigar_metrics[proportions[i]].clip(lower=0)

    read_name_to_cigar_metrics['best_hap'] = [get_best_hap_for_read(read_id, read_name_to_cigar_metrics) for read_id in
                                              read_name_to_cigar_metrics.index]
    return read_name_to_cigar_metrics


def get_best_hap_for_read(read_id, read_name_to_cigar_metrics):
    if read_id in read_name_to_cigar_metrics.index:
        hap1_prop = read_name_to_cigar_metrics.loc[read_id, 'hap1_prop']
        hap2_prop = read_name_to_cigar_metrics.loc[read_id, 'hap2_prop']
        if hap1_prop > hap2_prop:
            return 1
        else:
            return 2
    else:
        return 0


def cluster_haplotypes(gr1, gr2, read_name_to_cigar_metrics, output_file):
    print('Generating haplotype histogram metrics')

    print(read_name_to_cigar_metrics)

    legend_labels = ['Assigned to haplotype 1', 'Assigned to haplotype 2']
    legend_colors = ['gold', 'crimson']

    read_name_to_cigar_metrics_to_plot = read_name_to_cigar_metrics[['hap1_prop', 'hap2_prop']]

    row_colors = [legend_colors[0] if best_hap == 1 else legend_colors[1] if best_hap == 2 else 'black' for best_hap in
                  read_name_to_cigar_metrics['best_hap']]
    fig = plot_histogram_metrics(legend_colors, legend_labels, read_name_to_cigar_metrics_to_plot, row_colors,
                                'Read haplotype match metrics')

    print("saving figure to " + output_file)
    fig.set_size_inches(12, 10)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')


def plot_histogram_metrics(legend_colors, legend_labels, read_name_to_cigar_metrics_to_plot, row_colors, title):
    """
    Create histogram plots showing the distribution of haplotype metrics.
    This provides a sense of base-rate distance background against distances of events.
    """
    # Create subplots for histograms
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title, fontsize=16)
    
    # Get the data
    hap1_data = read_name_to_cigar_metrics_to_plot['hap1_prop']
    hap2_data = read_name_to_cigar_metrics_to_plot['hap2_prop']
    
    # Calculate the difference between haplotype proportions (distance metric)
    hap_diff = hap2_data - hap1_data
    
    # Separate data by assignment for colored histograms
    hap1_assigned = hap1_data[pd.Series(row_colors) == legend_colors[0]]
    hap1_assigned_hap2 = hap2_data[pd.Series(row_colors) == legend_colors[0]]
    hap2_assigned = hap1_data[pd.Series(row_colors) == legend_colors[1]]
    hap2_assigned_hap2 = hap2_data[pd.Series(row_colors) == legend_colors[1]]
    
    # Plot 1: Haplotype 1 proportions distribution
    axes[0, 0].hist(hap1_assigned, bins=30, alpha=0.7, color=legend_colors[0], label=legend_labels[0])
    axes[0, 0].hist(hap2_assigned, bins=30, alpha=0.7, color=legend_colors[1], label=legend_labels[1])
    axes[0, 0].set_xlabel('Haplotype 1 Proportion')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Haplotype 1 Proportions')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Haplotype 2 proportions distribution
    axes[0, 1].hist(hap1_assigned_hap2, bins=30, alpha=0.7, color=legend_colors[0], label=legend_labels[0])
    axes[0, 1].hist(hap2_assigned_hap2, bins=30, alpha=0.7, color=legend_colors[1], label=legend_labels[1])
    axes[0, 1].set_xlabel('Haplotype 2 Proportion')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of Haplotype 2 Proportions')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Difference between haplotype proportions (distance background)
    axes[1, 0].hist(hap_diff, bins=30, alpha=0.7, color='steelblue', edgecolor='black')
    axes[1, 0].set_xlabel('Haplotype Proportion Difference (H2 - H1)')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Haplotype Proportion Differences')
    axes[1, 0].axvline(x=0, color='red', linestyle='--', alpha=0.7, label='Equal proportions')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Combined scatter with density background
    axes[1, 1].scatter(hap1_data, hap2_data, c=row_colors, alpha=0.6, s=20)
    axes[1, 1].set_xlabel('Haplotype 1 Proportion')
    axes[1, 1].set_ylabel('Haplotype 2 Proportion')
    axes[1, 1].set_title('Haplotype Proportions Scatter Plot')
    
    # Add diagonal line for equal proportions
    max_val = max(hap1_data.max(), hap2_data.max())
    axes[1, 1].plot([0, max_val], [0, max_val], 'r--', alpha=0.7, label='Equal proportions')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    # Create custom legend for colors
    legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                                markersize=8, label=label) 
                     for color, label in zip(legend_colors, legend_labels)]
    fig.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(0.98, 0.95))
    
    plt.tight_layout()
    return fig
