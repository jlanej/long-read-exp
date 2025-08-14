#!/usr/bin/env python3
"""
Extract sequences from primary alignments starting at the query base that maps to a given
reference start coordinate, and keep reads that provide at least --length bases after that
position (insertions are preserved). Extraction goes from the anchor base to the end of the
read (including soft-clipped trailing bases).

Assumptions:
 - --start must correspond to a query position in the primary alignment (if start falls in a deletion
   there is no corresponding query base and the read is skipped).
 - --length refers to the minimum number of query bases extracted (including inserted bases).
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
    or None if no such query position exists in this alignment (e.g. start falls in a deletion).
    """
    try:
        pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
    except ValueError:
        return None

    for qpos, rpos in pairs:
        if rpos is None:
            continue
        if rpos == start_0:
            return qpos  # If start is in a deletion this will be None -> handled by caller
    return None

def collect_sequence_from_qstart_to_end(read, q_start):
    """
    Return the query sequence from q_start to the end of the read (includes soft-clipped bases).
    If read.query_sequence is None or q_start is out of bounds, return an empty string.
    """
    seq = read.query_sequence
    if seq is None:
        return ""

    if q_start is None:
        return ""

    # q_start should be an int; make sure it's within the query sequence bounds
    if q_start < 0 or q_start >= len(seq):
        return ""

    # Return from q_start to the end of the read (includes soft-clipped bases).
    return seq[q_start:]

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

        seq = collect_sequence_from_qstart_to_end(read, qstart)
        if not seq:
            continue
        if len(seq) < required_len:
            continue  # skip reads that don't provide at least required_len bases starting at the anchor

        rid = f"{sample_name}_{read.query_name}" if sample_name else read.query_name
        print(f">{rid}_{len(seq)}\n{seq}")

    bam.close()

if __name__ == "__main__":
    main()
