import argparse
import os

import numpy as np
import matplotlib.pyplot as plt


# adapted from https://medium.com/@anoopjohny2000/visualizing-sequence-similarity-with-dotplots-in-python-f5cf0ac8559f#:~:text=To%20create%20a%20Dot%20plot,with%20ungapped%20alignment%20and%20wordsize.
def dotplot(seqA, seqB, w, s):
    # Initialize the dotplot matrix with zeros
    dp = np.zeros((len(seqA), len(seqB)), dtype=int)
    window_b_cache = {}
    # Iterate over all positions in seqA and seqB
    for i in range(len(seqA)):
        if i % 100 == 0:
            print(f"processing base {i} of {len(seqA)}")

        # Compute the window around position i in seqA
        windowA = seqA[max(i - w, 0):min(i + w + 1, len(seqA))]
        for j in range(len(seqB)):
            # Compute the window around position j in seqB
            if j not in window_b_cache:
                window_b_cache[j] = seqB[max(j - w, 0):min(j + w + 1, len(seqB))]
            windowB = window_b_cache[j]
            # Count the number of matching symbols in the window
            matches = sum([1 for x, y in zip(windowA, windowB) if x == y])

            # Set the matrix element to 1 if the number of matches is at least s
            if matches >= s:
                dp[i, j] = 1
                # print(f"match at {i}, {j} with {matches} matches")

    return dp


def root_file_name_sans_dir(file_name):
    return file_name.split("/")[-1].split(".")[0]


def dotplot2Graphics(dp, labelA, labelB, heading, filename):
    # create a new figure
    fig, ax = plt.subplots()

    # plot the dots using a scatter plot
    rows, cols = np.where(dp)
    ax.scatter(cols, rows, marker='.', color='black')

    # set the labels and title
    ax.set_xlabel(root_file_name_sans_dir(labelB))
    ax.set_ylabel(root_file_name_sans_dir(labelA))
    ax.set_title(heading)

    # save the figure to a file and display it on screen
    plt.gcf().set_size_inches(25, 25)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    # plt.show()


def main():
    # Define the command line arguments
    parser = argparse.ArgumentParser(description="Generate a dotplot from two sequences in FASTA format.")
    parser.add_argument("w", type=int, help="the window size")
    parser.add_argument("s", type=int, help="the stringency")
    parser.add_argument("seqA", help="the filename of the first sequence in FASTA format")
    parser.add_argument("seqB", help="the filename of the second sequence in FASTA format")
    parser.add_argument("title", help="the title of the dotplot")
    parser.add_argument("output", help="the output filename for the dotplot")
    # argument for whether we should take the reverse complement of the first sequence
    parser.add_argument("--reverse", action="store_true", help="take the reverse complement of the first sequence")

    # Parse the command line arguments
    args = parser.parse_args()

    # Read the sequences from the input files
    with open(args.seqA) as fileA, open(args.seqB) as fileB:
        seqA = "".join([line.strip() for line in fileA if not line.startswith(">")])
        seqB = "".join([line.strip() for line in fileB if not line.startswith(">")])
    # make sure the sequences are upper case
    seqA = seqA.upper()
    seqB = seqB.upper()
    print("length of seqA: ", len(seqA))
    print("length of seqB: ", len(seqB))

    # Take the reverse complement of the first sequence if requested
    if args.reverse:
        seqA = seqA[::-1].translate(str.maketrans("ACGT", "TGCA"))
        print("T aking the reverse complement of the first sequence")
    cache_file = args.output + ".window." + str(args.w) + ".stringency." + str(args.s) + ".npy"
    # if the cache file does not exist ,gernerate the dp and save it

    if not os.path.exists(cache_file):
        dp = dotplot(seqA, seqB, args.w, args.s)
        np.save(cache_file, dp)
    else:
        dp = np.load(cache_file)
    # dp = dotplot(seqA, seqB, args.w, args.s)
    # # Generate the graphical dotplot
    dotplot2Graphics(dp, args.seqA, args.seqB, args.title, args.output + ".png")


if __name__ == "__main__":
    main()
