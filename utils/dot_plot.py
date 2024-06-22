import argparse
import numpy as np
import matplotlib.pyplot as plt
import pysam
import cigar
import wotplot


# https://github.com/fedarko/wotplot/blob/a61953c504bea7167cc06a29529d269fec09a9c6/wotplot/_matrix.py#L73C1-L77C60
#
#             -  2: k1 == k2, and ReverseComplement(k1) == k2
#             -  1: k1 == k2, and ReverseComplement(k1) != k2
#             - -1: k1 != k2, and ReverseComplement(k1) == k2
#             -  0: k1 != k2, and ReverseComplement(k1) != k2

COLORS = {
    2: "black",
    1: "green",
    -1: "red",
    0: "white"
}

def get_color_for_value(value):
    return COLORS[value]


def dotplot(seq1, seq2, w):
    return wotplot.DotPlotMatrix(seq1, seq2, w, binary=False, yorder="TB")


def root_file_name_sans_dir(file_name):
    return file_name.split("/")[-1].split(".")[0]


def get_start_index_from_label(label):
    #     if the first _ delimited entry does not start with chr, return 0
    if not label.split("_")[0].startswith("chr"):
        return 0
    return int(label.split("_")[1])


def get_sparse_subset_by_value(matrix, value):
    return matrix.mat == value

def dotplot2Graphics(dpo, label_x, label_y, heading, filename):
    # create a new figure
    dp = dpo.mat
    fig, ax = plt.subplots()

    # ax.spy(dp, marker='.', markersize=2, color='black', origin='lower')
    # loop over the possible values in the matrix and plot them
    for i in COLORS.keys():
        subset = get_sparse_subset_by_value(dpo, i)
        ax.spy(subset, marker='.', markersize=2, color=get_color_for_value(i), origin='lower')

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

    # ax.add_patch(plt.Rectangle((18000, 0), 100, dp.shape[0], fill=True, color='green', alpha=0.2))
    # ax.add_patch(plt.Rectangle((0, 18000), dp.shape[1], 100, fill=True, color='blue', alpha=0.2))
    # save the figure to a file and display it on screen
    plt.gcf().set_size_inches(25, 25)
    print("saving png file: ", filename)
    plt.savefig(filename, dpi=300, bbox_inches='tight')


def get_png_file_for_read(read, k, output):
    sanitized_read_name = read.query_name.replace("/", "_")
    output_read = output + "." + sanitized_read_name
    cache_file = output_read + ".k." + str(k) + ".png"
    return cache_file


def has_hard_clipping(read):
    #     if the read has hard clipping
    cr = cigar.Cigar(read.cigarstring)
    for c in cr.items():
        if c[1] == "H":
            return True
    return False


def non_primary_alignment(read):
    return read.is_secondary or read.is_supplementary


def use_read(read):
    return not has_hard_clipping(read) and not non_primary_alignment(read)


def dot_fasta_vs_fasta(reference_seq_file, compSeq, k, output):
    with open(reference_seq_file) as fileA:
        seq_ref = "".join([line.strip() for line in fileA if not line.startswith(">")])
    seq_ref = seq_ref.upper()
    with open(compSeq) as fileB:
        seq_comp = "".join([line.strip() for line in fileB if not line.startswith(">")])
    seq_comp = seq_comp.upper()

    print("length of ref seq: ", len(seq_ref))
    print("length of comp seq: ", len(seq_comp))
    dp = dotplot(seq_ref, seq_comp, k)
    dotplot2Graphics(dp, reference_seq_file, compSeq, root_file_name_sans_dir(reference_seq_file) + " vs " + root_file_name_sans_dir(compSeq),
                     output + ".k." + str(k) + ".png")


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


def main():
    # Define the command line arguments
    parser = argparse.ArgumentParser(description="Generate a dotplot from two sequences in FASTA format.")
    parser.add_argument("--k", type=int, help="kmer size")
    parser.add_argument("--bam", help="the bam file of sequences")
    parser.add_argument("--reference_seq", help="the filename of reference in FASTA format")
    parser.add_argument("--compSeq", help="the filename of a second sequence in FASTA format")
    parser.add_argument("--output", help="the root output filename for the dotplot")
    args = parser.parse_args()

    # if bam and compSeq arguments are both provided, exit
    if args.bam and args.compSeq:
        print("Please provide either a BAM file or a second sequence (compSeq) in FASTA format, not both.")
        return
    # if compSeq argument is provided, generate a dotplot from two sequences in FASTA format
    if args.compSeq:
        dot_fasta_vs_fasta(args.reference_seq, args.compSeq, args.k, args.output)
        return

    a = pysam.AlignmentFile(args.bam, "rb")
    for read in a:
        if not use_read(read):
            print("skipping read ", read.query_name, " because it is not primary alignment or has hard clipping")
            continue
        dp = dot_read(args.reference_seq, read, args.k)
        dotplot2Graphics(dp, args.reference_seq, read.query_name,
                         read.query_name + " vs " + root_file_name_sans_dir(args.reference_seq),
                         get_png_file_for_read(read, args.k, args.output))


if __name__ == "__main__":
    main()
