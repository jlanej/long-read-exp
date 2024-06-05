import argparse
import numpy as np
import matplotlib.pyplot as plt
# adapted from https://medium.com/@anoopjohny2000/visualizing-sequence-similarity-with-dotplots-in-python-f5cf0ac8559f#:~:text=To%20create%20a%20Dot%20plot,with%20ungapped%20alignment%20and%20wordsize.
def dotplot(seqA, seqB, w, s):
    # Initialize the dotplot matrix with zeros
    dp = np.zeros((len(seqA), len(seqB)), dtype=int)

    # Iterate over all positions in seqA and seqB
    for i in range(len(seqA)):
        if i % 100 == 0:
            print(f"processing base {i} of {len(seqA)}")

        # Compute the window around position i in seqA
        windowA = seqA[max(i - w, 0):min(i + w + 1, len(seqA))]
        for j in range(len(seqB)):
            # Compute the window around position j in seqB
            windowB = seqB[max(j-w, 0):min(j+w+1, len(seqB))]
            # Count the number of matching symbols in the window
            matches = sum([1 for x, y in zip(windowA, windowB) if x == y])

            # Set the matrix element to 1 if the number of matches is at least s
            if matches >= s:
                dp[i,j] = 1
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

    # set the tick positions and labels
    # xticks = np.arange(0.5, dp.shape[1], 1)
    # yticks = np.arange(0.5, dp.shape[0], 1)
    # ax.set_xticks(xticks)
    # ax.set_yticks(yticks)
    # ax.set_xticklabels(np.arange(1, dp.shape[1] + 1))
    # ax.set_yticklabels(np.arange(1, dp.shape[0] + 1)[::-1])

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

    # Parse the command line arguments
    args = parser.parse_args()

    # Read the sequences from the input files
    with open(args.seqA) as fileA, open(args.seqB) as fileB:
        seqA = "".join([line.strip() for line in fileA if not line.startswith(">")])
        seqB = "".join([line.strip() for line in fileB if not line.startswith(">")])
    print(len(seqA), len(seqB))
    print("length of seq from file",args.seqA, " is ", len(seqA))
    print("length of seq from file",args.seqB, " is ", len(seqB))
    # trim seqb to 1000bases
    # seqB = seqB[:1000]
    # Generate the dotplot matrix
    print("creating dotplot")
    dp = dotplot(seqA, seqB, args.w, args.s)
    # Generate the graphical dotplot
    dotplot2Graphics(dp, args.seqA, args.seqB, args.title, args.output + ".png")

if __name__ == "__main__":
    main()