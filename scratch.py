

def harmonize_cigars(gr,gr_reduce):
    all_start = gr_reduce.get_start()[0]
    all_end = gr_reduce.get_end()[0]
    gr.mcols.set_column('original_cigar', gr.mcols.get_column('cigar'), in_place=True)
    for i in range(len(gr)):
        cigars = cigar.Cigar(gr.mcols.get_column('cigar')[i])
        start = gr.get_start()[i]
        num_soft_clipping = start - all_start
        effective_end = start + cigars.reference_length()
        if start > all_start:
            cigars = cigar.Cigar(str(num_soft_clipping) + 'S'+str(cigars))
        if effective_end < all_end:
            print(cigars)
            print(len(cigars))
            cigars = cigar.Cigar(str(cigars) + str(all_end - effective_end) + 'S')
            print(cigars)
            print(len(cigars))

        gr.mcols.get_column('cigar')[i] = str(cigars)
        if len(cigars) != all_end - all_start:
            sys.stderr.write('Error: Cigar length does not match reference genome\n')
            sys.stderr.write(str(len(cigars)) + ' ' + str(all_end - all_start) + '\n')
            sys.exit(1)
    # make sure the cigar strings are the same length

    return gr

# makes all read sequences the same length by adding Ns to the beginning or end of the sequence

def harmonize_sequences(gr, gr_reduce):
    all_start = gr_reduce.get_start()[0]
    all_end = gr_reduce.get_end()[0]
    gr.mcols.set_column('original_sequence', gr.mcols.get_column('sequence'), in_place=True)

    for i in range(len(gr)):
        cigars = cigar.Cigar(gr.mcols.get_column('cigar')[i])
        start = gr.get_start()[i]
        current_sequence = gr.mcols.get_column('sequence')[i]
        num_Ns = start - all_start
        if start > all_start:
            current_sequence = '.' * num_Ns + current_sequence

        indexToInsert = num_Ns
        for c in cigars.items():
            #     if the operation consumes the reference genome but not the read, add "." to the read sequence
            #     to match the reference genome
            if c[1] in ['D']:
                current_sequence = current_sequence[:indexToInsert] + '.' * c[0] + current_sequence[indexToInsert:]
                indexToInsert += c[0]
        effective_end = len(current_sequence)

        if effective_end < all_end:
            current_sequence = current_sequence + '.' * (all_end - effective_end)
        if len(current_sequence) < all_end - all_start:
            sys.stderr.write('Error: Sequence length too short\n')
            # print the sequence length and the expected length
            print(effective_end)
            print(len(current_sequence), all_end - all_start, file=sys.stderr)
            sys.exit(1)
        # trim the sequence to 10000 bases on either side of the middle of the full range
        middle = (all_end - all_start) // 2
        # current_sequence = current_sequence[middle - 10000:middle + 10000]
        gr.mcols.get_column('sequence')[i] = current_sequence

    return gr


def create_harmonized_sequence(alignment_file):
    reads = extract_reads(alignment_file)
    gr, gr_reduce = convert_to_range(reads)
    # gr = harmonize_sequences(gr, gr_reduce)
    harmonize_cigars(gr, gr_reduce)
    return gr



# create a plot of the harmonized sequences to check that they are aligned
def plot_harmonized_sequence(gr):
    fig, ax = plt.subplots()
    for i in range(len(gr)):
        ax.text(0, i, gr.mcols.get_column('read_name')[i], fontsize=12)
        ax.text(1, i, gr.mcols.get_column('sequence')[i], fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(gr))
    plt.show()
