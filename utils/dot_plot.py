import argparse
import os
import sys

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pysam
import gzip
import cigar
import wotplot

from multiprocessing.dummy import Pool as ThreadPool


# adapted from https://medium.com/@anoopjohny2000/visualizing-sequence-similarity-with-dotplots-in-python-f5cf0ac8559f#:~:text=To%20create%20a%20Dot%20plot,with%20ungapped%20alignment%20and%20wordsize.
# def dotplot(seqA, seqB, w):
#     # Initialize the dotplot matrix with zeros
#     dp = np.zeros((len(seqA), len(seqB)), dtype=int)
#     window_b_cache = {}
#     # Iterate over all positions in seqA and seqB
#     for i in range(len(seqA)):
#         if i % 100 == 0:
#             print(f"processing base {i} of {len(seqA)}")
#
#         # Compute the window around position i in seqA
#         windowA = seqA[max(i - w, 0):min(i + w + 1, len(seqA))]
#         for j in range(len(seqB)):
#             # Compute the window around position j in seqB
#             if j not in window_b_cache:
#                 window_b_cache[j] = seqB[max(j - w, 0):min(j + w + 1, len(seqB))]
#             windowB = window_b_cache[j]
#             # Count the number of matching symbols in the window
#             matches = sum([1 for x, y in zip(windowA, windowB) if x == y])
#
#             # Set the matrix element to 1 if the number of matches is at least s
#             dp[i, j] = matches
#
#             # print(f"match at {i}, {j} with {matches} matches")
#
#     return dp
def dotplot(seqA, seqB, w):
    return wotplot.DotPlotMatrix(seqA, seqB, w, binary=False).mat

def wot_dotplot(seqA, seqB,k):
    # m = wotplot.DotPlotMatrix(s1, s2, k, binary=False)
    return wotplot.DotPlotMatrix(seqA, seqB, k, binary=False).mat

def root_file_name_sans_dir(file_name):
    return file_name.split("/")[-1].split(".")[0]


def get_start_index_from_label(label):
#     if the first _ delimited entry does not start with chr, return 0
    if not label.split("_")[0].startswith("chr"):
        return 0
    return int(label.split("_")[1])
def dotplot2Graphics(dp, s, labelA, labelB, heading, filename):
    # create a new figure
    fig, ax = plt.subplots()

    # plot the dots using a scatter plot
    # need where to get the indices of the entries in dp that greater than or equal to s
    rows, cols = np.where(dp >= s)
    # rows, cols = np.where(dp)
    ax.scatter(cols, rows, marker='.', color='black')

    labelAUse=root_file_name_sans_dir(labelA)
    labelBUse=root_file_name_sans_dir(labelB)


    # determine labelEvery to have 10 labels on the x axis
    labelEvery=dp.shape[1]//10
    ax.set_xticks(np.arange(0, dp.shape[1], labelEvery))
    startx=get_start_index_from_label(labelBUse)
    ax.set_xticklabels(np.arange(startx, startx + dp.shape[1], labelEvery))

    # determine labelEvery to have 10 labels on the y axis
    labelEvery=dp.shape[0]//10
    ax.set_yticks(np.arange(0, dp.shape[0], labelEvery))
    starty=get_start_index_from_label(labelAUse)
    ax.set_yticklabels(np.arange(starty, starty + dp.shape[0], labelEvery))

    # set the labels and title
    ax.set_xlabel(root_file_name_sans_dir(labelBUse))
    ax.set_ylabel(root_file_name_sans_dir(labelAUse))
    ax.set_title(heading)

    # add a transparent rectangle to highlight the region of interest from x=181 to x=13993 and the max y index, colored green
    ax.add_patch(plt.Rectangle((18000, 0),100, dp.shape[0], fill=True, color='green', alpha=0.2))
    ax.add_patch(plt.Rectangle((0, 18000),dp.shape[1], 100, fill=True, color='blue', alpha=0.2))

    # ax.add_patch(plt.Rectangle((15877, 0), 17586 - 15877, dp.shape[0], fill=True, color='purple', alpha=0.2))
    # ax.add_patch(plt.Rectangle((13994, 0), 17586 - 15877, dp.shape[0], fill=True, color='yellow', alpha=0.2))
    # ax.add_patch(plt.Rectangle((0, 0), dp.shape[1], 13812, fill=True, color='green', alpha=0.2))
    # ax.add_patch(plt.Rectangle((17980, 0), 10, 1710, fill=True, color='blue', alpha=0.2))

    # save the figure to a file and display it on screen
    plt.gcf().set_size_inches(25, 25)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    # plt.show()


def get_cache_file_for_read(read, w, output):
    sanitized_read_name = read.query_name.replace("/", "_")
    output_read = output + "." + sanitized_read_name
    cache_file = output_read + ".window." + str(w) + ".npy.gz"
    return cache_file


def get_png_file_for_read(read, w, s, output):
    sanitized_read_name = read.query_name.replace("/", "_")
    output_read = output + "." + sanitized_read_name
    cache_file = output_read + ".window." + str(w) + ".s." + str(s) + ".png"
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



def dot_fasta_vs_fasta(fastaA, fastaB, w, s, output):
    with open(fastaA) as fileA:
        seqA = "".join([line.strip() for line in fileA if not line.startswith(">")])
    seqA = seqA.upper()
    with open(fastaB) as fileB:
        seqB = "".join([line.strip() for line in fileB if not line.startswith(">")])
    seqB = seqB.upper()

    print("length of seqA: ", len(seqA))
    print("length of seqB: ", len(seqB))

    # if the cache file does not exist ,generate the dp and save it
    cache_file = output + ".window." + str(w)
    if not os.path.exists(cache_file+".npz"):
        dp = dotplot(seqA, seqB, w)
        # f = gzip.GzipFile(cache_file, "w")
        # np.save(f, dp)
        sp.sparse.save_npz(cache_file, dp)
        # f.close()
    else:
        print(cache_file, " already exists")
    # dp = load_dot(cache_file)
    dp = sp.sparse.load_npz(cache_file+".npz")
    dotplot2Graphics(dp, s, fastaA, fastaB, root_file_name_sans_dir(fastaA) + " vs " + root_file_name_sans_dir(fastaB), output + ".window." + str(w) + ".s." + str(s) + ".png")

def prep_dot(read, reference_seq_file, w, output):
    if not use_read(read):
        print("skipping read ", read.query_name, " because it is not primary alignment or has hard clipping")
        return None
    cache_file = get_cache_file_for_read(read, w, output)
    with open(reference_seq_file) as fileB:
        reference_seq = "".join([line.strip() for line in fileB if not line.startswith(">")])
    reference_seq = reference_seq.upper()

    seqA = read.seq

    print("read name: ", read.query_name)
    print("cache file: ", cache_file)
    print("length of seqA: ", len(seqA))
    print("length of seqB: ", len(reference_seq))
    # Take the reverse complement of the first sequence if requested
    if read.is_reverse:
        # seqA = seqA[::-1].translate(str.maketrans("ACGT", "TGCA"))
        print("Taking the reverse complement of the first sequence")

    # if the cache file does not exist ,generate the dp and save it
    if not os.path.exists(cache_file):
        dp = dotplot(seqA, reference_seq, w)
        f = gzip.GzipFile(cache_file, "w")
        np.save(f, dp)

        f.close()
    else:
        print(cache_file, " already exists")
    return cache_file


def load_dot(cache_file):
    print("loading ", cache_file)
    f = gzip.GzipFile(cache_file, "r")
    dp = np.load(f)
    f.close()
    return dp

def main():
    # Define the command line arguments
    parser = argparse.ArgumentParser(description="Generate a dotplot from two sequences in FASTA format.")
    parser.add_argument("--w", type=int, help="the window size")
    parser.add_argument("--s", type=int, help="the stringency")
    parser.add_argument("--bam", help="the bam file of sequences")
    parser.add_argument("--reference_seq", help="the filename of reference in FASTA format")
    parser.add_argument("--compSeq", help="the filename of a second sequence in FASTA format")
    parser.add_argument("--output", help="the root output filename for the dotplot")
    parser.add_argument("--threads", type=int, help="the number of threads to use", default=8)
    args = parser.parse_args()

    # if compSeq argument is provided, generate a dotplot from two sequences in FASTA format
    if args.compSeq:
        dot_fasta_vs_fasta(args.reference_seq, args.compSeq, args.w, args.s, args.output)
        return

    threads = args.threads
    print("using ", threads, " threads")
    pool = ThreadPool(threads)
    pool.map(lambda r: prep_dot(r, args.reference_seq, args.w, args.output), pysam.AlignmentFile(args.bam, "rb"))
    pool.close()
    a = pysam.AlignmentFile(args.bam, "rb")
    for read in a:
        if not use_read(read):
            print("skipping read ", read.query_name, " because it is not primary alignment or has hard clipping")
            continue
        cache_file = get_cache_file_for_read(read, args.w, args.output)
        dp = load_dot(cache_file)
        dotplot2Graphics(dp, args.s, read.query_name, args.reference_seq, read.query_name + " vs " +root_file_name_sans_dir(args.reference_seq),
                         get_png_file_for_read(read, args.w, args.s, args.output))


if __name__ == "__main__":
    main()
