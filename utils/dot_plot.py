import argparse
import numpy as np
import matplotlib.pyplot as plt
import pysam
import cigar
import wotplot
from matplotlib.legend import Legend
import lr_utils
# https://github.com/fedarko/wotplot/blob/a61953c504bea7167cc06a29529d269fec09a9c6/wotplot/_matrix.py#L73C1-L77C60
#
#             -  2: k1 == k2, and ReverseComplement(k1) == k2
#             -  1: k1 == k2, and ReverseComplement(k1) != k2
#             - -1: k1 != k2, and ReverseComplement(k1) == k2
#             -  0: k1 != k2, and ReverseComplement(k1) != k2
DOT_COLORS = {
    2: ["black", "palindrome"],
    1: ["green", "forward match"],
    -1: ["red", "reverse complement match"],
    0: ["white", "no match"]
}

CIGAR_COLORS = {
    "M": "green",
    "I": "blue",
    "D": "red",
    "S": "gray",
    "H": "gray",
    "P": "gray",
    "=": "green",
    "X": "green"
}


def dotplot(seq1, seq2, w):
    return wotplot.DotPlotMatrix(seq1, seq2, w, binary=False, yorder="TB", verbose=True)


def root_file_name_sans_dir(file_name):
    return file_name.split("/")[-1].split(".")[0]


def get_start_index_from_label(label):
    #     if the first _ delimited entry does not start with chr, return 0
    if not label.split("_")[0].startswith("chr"):
        return 0
    return int(label.split("_")[1])


def get_chr_from_label(label):
    if not label.split("_")[0].startswith("chr"):
        return None
    return label.split("_")[0]


def get_ref_loc_from_label(label):
    if not label.split("_")[0].startswith("chr"):
        return [None, None, None]
    return [label.split("_")[0], int(label.split("_")[1]), int(label.split("_")[2])]


def get_sparse_subset_by_value(matrix, value):
    return matrix.mat == value


def plot_dot(dpo, ax, label_x, label_y, heading, marker_size):
    # create a new figure
    dp = dpo.mat
    # loop over the possible values in the matrix and plot them
    for i in DOT_COLORS.keys():
        # skip 0
        if i == 0:
            continue
        subset = get_sparse_subset_by_value(dpo, i)
        ax.spy(subset, marker='.', markersize=marker_size, color=DOT_COLORS[i][0], origin='lower', alpha=0.5)

    # set the x and y-axis labels
    label_x_use = root_file_name_sans_dir(label_x)
    label_y_use = root_file_name_sans_dir(label_y)

    # determine labelEvery to have 10 labels on the x axis
    label_every = dp.shape[1] // 10
    ax.set_xticks(np.arange(0, dp.shape[1], label_every))
    start_x = get_start_index_from_label(label_x_use)
    ax.set_xticklabels(np.arange(start_x, start_x + dp.shape[1], label_every))
    #
    # # determine labelEvery to have 10 labels on the y axis
    label_every = dp.shape[0] // 10
    ax.set_yticks(np.arange(0, dp.shape[0], label_every))
    start_y = get_start_index_from_label(label_y_use)
    ax.set_yticklabels(np.arange(start_y, start_y + dp.shape[0], label_every))

    # set the labels and title
    ax.set_xlabel(root_file_name_sans_dir(label_x_use))
    ax.set_ylabel(root_file_name_sans_dir(label_y_use))
    ax.set_title(heading)
    return ax


def save_plot(ax, filename):
    plt.gcf().set_size_inches(25, 25)
    print("saving png file: ", filename)
    plt.savefig(filename, dpi=300, bbox_inches='tight')


def get_png_file_for_read(read, k, output, cigar=False):
    sanitized_read_name = read.query_name.replace("/", "_")
    output_read = output + "." + sanitized_read_name
    cache_file = output_read + ".k." + str(k) + (".cigar" if cigar else "") + ".png"
    return cache_file


def has_hard_clipping(read):
    #     if the read has hard clipping
    cr = cigar.Cigar(read.cigarstring)
    for c in cr.items():
        if c[1] == "H":
            return True
    return False


def use_read(read):
    return not has_hard_clipping(read)


def dot_fasta_vs_fasta(reference_seq_file, compSeq, k, output, marker_size, create_legend_plot=False):
    with open(reference_seq_file) as fileA:
        seq_ref = "".join([line.strip() for line in fileA if not line.startswith(">")])
    seq_ref = seq_ref.upper()
    with open(compSeq) as fileB:
        seq_comp = "".join([line.strip() for line in fileB if not line.startswith(">")])
    seq_comp = seq_comp.upper()

    print("length of ref seq: ", len(seq_ref))
    print("length of comp seq: ", len(seq_comp))
    dp = dotplot(seq_ref, seq_comp, k)
    fig, ax = plt.subplots()

    ax = plot_dot(dp, ax, reference_seq_file, compSeq,
                  root_file_name_sans_dir(reference_seq_file) + " vs " + root_file_name_sans_dir(
                      compSeq) + "\nk=" + str(k), marker_size)
    save_plot(ax, output + ".k." + str(k) + ".png")
    if create_legend_plot:
        # Add a legend describing the colors
        add_legend(ax)
        save_plot(ax, output + ".k." + str(k) + ".legend.png")


def add_legend(ax):
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label=DOT_COLORS[i][1], markerfacecolor=DOT_COLORS[i][0],
                   markersize=10) for i in DOT_COLORS.keys() if i != 0]
    ax.legend(handles=legend_elements, loc='upper left', title="Dot colors", bbox_to_anchor=(0.25, -0.1))


def dot_read(reference_seq_file, read, k):
    if not use_read(read):
        print("skipping read ", read.query_name, " because it is not primary alignment or has hard clipping")
        return None
    with open(reference_seq_file) as fileB:
        reference_seq = "".join([line.strip() for line in fileB if not line.startswith(">")])
    reference_seq = reference_seq.upper()

    print("read name: ", read.query_name)
    print("length of read: ", len(read.seq))
    print("length of ref: ", len(reference_seq))

    return dotplot(reference_seq, read.seq, k)


def get_base_alignment(read):
    cr = cigar.Cigar(read.cigarstring)
    return [read.reference_name, read.reference_start, read.reference_end, cr]


def parse_SA_tag(tag):
    sa_all = tag.split(";")
    sas = []
    for sa in sa_all:
        if sa == "":
            continue
        split_sa = sa.split(",")
        chr = split_sa[0]
        start = int(split_sa[1])
        cigarS = cigar.Cigar(split_sa[3])
        end = start + cigarS.reference_length()
        sas.append([chr, start, end, cigarS])
    return sas


def get_all_alignments(read):
    # store a string of chr:start-end for each alignment and the cigar string associated with it
    alignments = []
    alignments.append(get_base_alignment(read))
    for tag in read.tags:
        if tag[0] == "SA":
            sas = parse_SA_tag(tag[1])
            for sa in sas:
                alignments.append(sa)
    alignments = sort_alignments_by_left_clip(alignments)
    return alignments


def get_amount_of_left_soft_hard_clipping(alignment):
    left_clip = 0
    cr = alignment[3]
    for c in cr.items():
        if c[1] in ["H", "S"]:
            left_clip += c[0]
        else:
            break
    return left_clip


def sort_alignments_by_left_clip(alignments):
    return sorted(alignments, key=lambda x: get_amount_of_left_soft_hard_clipping(x))


# "N": "white",

def collapse_cigar_colors():
    #     group CIGAR_COLORS that are the same color
    color_groups = {}
    for c in CIGAR_COLORS.items():
        if c[1] not in color_groups:
            color_groups[c[1]] = ""
        else:
            color_groups[c[1]] += ","
        color_groups[c[1]] += c[0]
    return color_groups


def get_label_for_alignment(alignment, number):
    return alignment[0] + ":" + str(alignment[1]) + "-" + str(alignment[2]) + " (#" + str(number) + ")"


def add_cigar_to_fig(ax, read, min_indel, ref_loc):
    alignments = get_all_alignments(read)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    width = (xmax - xmin) * 0.1
    per_alignment_width = width / len(alignments)
    rect_x = xmin - width

    # set the new xmin
    ax.set_xlim(rect_x, xmax)
    # draw a vertical line at the old xmin
    ax.axvline(x=xmin, color='black', linestyle='solid')
    alignment_number = 0

    labels_to_add = []
    ticks_to_add = []
    for alignment in alignments:
        alignment_x = rect_x + per_alignment_width * alignment_number
        ax.axvline(x=alignment_x + per_alignment_width, color='black', linestyle='dashed')
        alignment_number += 1
        labels_to_add.append(get_label_for_alignment(alignment, alignment_number))
        ticks_to_add.append(alignment_x + per_alignment_width / 2)
        same_chr = ref_loc[0] == alignment[0]
        read_index = 0
        ref_index = 0
        if not ref_loc[1] ==None:
            ref_index = 0 - (ref_loc[1] - alignment[1])
        for c in alignment[3].items():
            if c[1] in ["M", "X", "="]:
                ax.add_patch(
                    plt.Rectangle((alignment_x, read_index), per_alignment_width, c[0], fill=True,
                                  color=CIGAR_COLORS[c[1]], alpha=0.25,
                                  edgecolor=None))
                read_index += c[0]
                ref_index += c[0]
            elif c[1] in ["H", "S"]:
                ax.add_patch(
                    plt.Rectangle((alignment_x, read_index), per_alignment_width, c[0], fill=True,
                                  color=CIGAR_COLORS[c[1]], alpha=0.25,
                                  edgecolor=None))
                read_index += c[0]
            elif c[1] in ["D"]:
                if c[0] > min_indel and same_chr:
                    ax.add_patch(
                        plt.Rectangle((ref_index, 0), c[0], ymax, fill=True, color=CIGAR_COLORS[c[1]],
                                      alpha=0.15))
                ref_index += c[0]
                # ax.add_patch(plt.Rectangle((rect_x, read_index), xmax, c[0], fill=True, color=CIGAR_COLORS[c[1]], alpha=0.25))
            elif c[1] in ["I"]:
                if c[0] > min_indel:
                    ax.add_patch(
                        plt.Rectangle((alignment_x, read_index), per_alignment_width, c[0], fill=True,
                                      color=CIGAR_COLORS[c[1]],
                                      alpha=0.25))
                    if same_chr:
                        ax.add_patch(
                            plt.Rectangle((xmin, read_index), xmax, c[0], fill=True,
                                          color=CIGAR_COLORS[c[1]],
                                          alpha=0.15))
                read_index += c[0]
            else:
                print("unknown cigar: ", c)

    current_labels = ax.get_xticklabels()
    current_ticks = ax.get_xticks()
    # add the new x ticks at the beginning of the list
    ax.set_xticks(np.append(ticks_to_add, current_ticks))
    ax.set_xticklabels(np.append(labels_to_add, current_labels))
    # rotate the labels, but only the new ones
    # ax.set_xticklabels(current_labels, rotation=45, ha='right')
    # loop over the ticks and rotate the new ones
    for tick in ax.get_xticklabels():
        if tick.get_text() in labels_to_add:
            tick.set_rotation(45)
            tick.set_ha('right')
        else:
            tick.set_rotation(-45)
            tick.set_ha('left')
    # plt.xticks(rotation=45, ha='right')
    return ax


def main():
    # Define the command line arguments
    parser = argparse.ArgumentParser(description="Generate a dotplot from two sequences in FASTA format.")
    parser.add_argument("--k", type=int, help="kmer size")
    parser.add_argument("--bam", help="the bam file of sequences")
    parser.add_argument("--reference_seq", help="the filename of reference in FASTA format")
    parser.add_argument("--compSeq", help="the filename of a second sequence in FASTA format")
    parser.add_argument("--output", help="the root output filename for the dotplot")
    parser.add_argument("--marker_size", type=float, default=1, help="the size of the markers in the dotplot")
    parser.add_argument("--min_indel", type=int, default=10, help="minimum indel size to show in dotplot")
    args = parser.parse_args()

    # if bam and compSeq arguments are both provided, exit
    if args.bam and args.compSeq:
        print("Please provide either a BAM file or a second sequence (compSeq) in FASTA format, not both.")
        return
    # if compSeq argument is provided, generate a dotplot from two sequences in FASTA format
    if args.compSeq:
        dot_fasta_vs_fasta(args.reference_seq, args.compSeq, args.k, args.output, args.marker_size)
        return
    elif args.bam:
        dot_fasta_vs_fasta(args.reference_seq, args.reference_seq, args.k, args.output+".ref_v_ref", args.marker_size)
        process_bam(args)
        return
    else:
        print("Please provide either a BAM file or a second sequence (compSeq) in FASTA format.")
        return


def process_bam(args,create_legend_plot=False):
    a = pysam.AlignmentFile(args.bam, "rb")
    for read in a:
        if has_hard_clipping(read):
            print("skipping read ", read.query_name, " because it has hard clipping")
            continue
        dp = dot_read(args.reference_seq, read, args.k)
        fig, ax = plt.subplots()
        fig.tight_layout()
        plt.rcParams.update({'font.size': 22})
        ax = plot_dot(dp, ax, args.reference_seq, "read sequence index",
                      read.query_name + " vs " + root_file_name_sans_dir(args.reference_seq) + "\nk=" + str(args.k),
                      args.marker_size)
        save_plot(ax, get_png_file_for_read(read, args.k, args.output))
        ax = add_cigar_to_fig(ax, read, args.min_indel,
                              get_ref_loc_from_label(root_file_name_sans_dir(args.reference_seq)))

        save_plot(ax, get_png_file_for_read(read, args.k, args.output, cigar=True))

        if create_legend_plot:
            add_legend(ax)
            # add a color legend for the CIGAR colors
            groups = collapse_cigar_colors()
            keys = groups.keys()
            legend_elements = [
                plt.Line2D([0], [0], marker='o', color='w', label=groups[i], markerfacecolor=i,
                           markersize=10) for i in keys]
            labels = [groups[i] for i in keys]
            leg = Legend(ax, legend_elements, loc='upper right', title="CIGAR colors", labels=labels,
                         bbox_to_anchor=(0.75, -0.1))
            ax.add_artist(leg)
            save_plot(ax, get_png_file_for_read(read, args.k, args.output, cigar=True).replace(".png", ".legend.png"))
        plt.close()


if __name__ == "__main__":
    main()
