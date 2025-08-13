#!/usr/bin/env python3
"""
Extract sequences from primary alignments starting at the query base that maps to a given
reference start coordinate, and keep reads that provide at least --length bases after that
position (insertions are preserved).

Assumptions:
 - --start must correspond to a query position in the primary alignment (if start falls in a deletion
   there is no corresponding query base and the read is skipped).
 - --length refers to the number of query bases extracted (including inserted bases).
 - Only primary alignments are used; supplementary and secondary alignments are skipped.
"""
import pysam
import argparse

def get_sample_name(bam):
    samples = {rg.get("SM") for rg in bam.header.get("RG", []) if "SM" in rg}
    samples.discard(None)
    if len(samples) == 1:
        return samples.pop()
    elif len(samples) > 1:
        raise ValueError(f"Multiple sample names in BAM header: {samples}")
    return None

def find_qstart_for_ref_start(read, start_0):
    """
    Return the query position (qpos) that aligns to reference coordinate start_0 (exact equality),
    or None if no such query position exists in this alignment.
    """
    try:
        pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
    except ValueError:
        return None

    for qpos, rpos in pairs:
        if rpos is None:
            continue
        if rpos == start_0:
            return qpos  # may still be None, but get_aligned_pairs gives qpos for matches/insertions
    return None

def collect_sequence_after_qstart(read, q_start, length):
    """
    Build the sequence in query order starting at q_start, using the read.query_sequence.
    We include any query positions >= q_start that are part of the alignment (including insertions).
    Stop once we collected >= length bases.
    """
    if read.query_sequence is None:
        return ""

    # Build a set/list of qpos that appear in the aligned pairs (so we avoid taking soft-clipped trailing bases)
    try:
        pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
    except ValueError:
        return ""

    # Collect all qpos that map in any way (including insertions: rpos can be None)
    qpos_set = set()
    for qpos, rpos in pairs:
        if qpos is not None:
            qpos_set.add(qpos)

    # Create sorted list of qpos in ascending query order that are >= q_start
    qpos_list = sorted([q for q in qpos_set if q >= q_start])

    # Build sequence until we reach desired length
    seq_chunks = []
    collected = 0
    L = len(read.query_sequence)
    for q in qpos_list:
        if q < 0 or q >= L:
            continue
        seq_chunks.append(read.query_sequence[q])
        collected += 1
        if collected >= length:
            break

    return "".join(seq_chunks)

def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences from primary alignments starting at the reference start coordinate."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM/CRAM (indexed)")
    parser.add_argument("-c", "--chrom", required=True, help="Reference contig name")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start coordinate (1-based)")
    parser.add_argument("-l", "--length", type=int, required=True, help="Minimum length of bases to extract from the read (query-bases)")
    parser.add_argument("--sample", help="Optional sample name to prepend to read IDs (overrides BAM header)")
    args = parser.parse_args()

    start_0 = args.start - 1
    required_len = args.length

    bam = pysam.AlignmentFile(args.bam, "rb")
    sample_name = args.sample if args.sample else get_sample_name(bam)

    # iterate reads that overlap the start position
    for read in bam.fetch(args.chrom, start_0, start_0 + 1):
        # skip unmapped or non-primary
        if read.is_unmapped or read.is_supplementary or read.is_secondary:
            continue

        # find query position that maps exactly to start_0
        qstart = find_qstart_for_ref_start(read, start_0)
        if qstart is None:
            # no query base aligns exactly to the reference start (start in deletion or not aligned here)
            continue

        seq = collect_sequence_after_qstart(read, qstart, required_len)
        if not seq:
            continue
        if len(seq) < required_len:
            continue  # skip reads that don't provide at least required_len bases starting at the anchor

        rid = f"{sample_name}_{read.query_name}" if sample_name else read.query_name
        print(f">{rid}_{len(seq)}\n{seq}")

    bam.close()

if __name__ == "__main__":
    main()
