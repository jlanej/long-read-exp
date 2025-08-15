#!/usr/bin/env python3
"""
Extract sequences from alignments (primary + supplemental) starting at the query base that maps
to a given reference start coordinate, and keep reads that provide at least --length bases after
that position (insertions are preserved). Extraction goes from the anchor base to the end of the
read (including soft-clipped trailing bases).

Assumptions:
 - --start must correspond to a query position in the alignment (if start falls in a deletion,
   there is no corresponding query base and the read is skipped).
 - --length refers to the minimum number of query bases extracted (including inserted bases).
 - Both primary and supplemental alignments are used; secondary alignments are skipped.
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
    """Return the query position that aligns exactly to reference coordinate start_0."""
    try:
        pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
    except ValueError:
        return None
    for qpos, rpos in pairs:
        if rpos == start_0:
            return qpos
    return None

def collect_sequence_from_qstart_to_end(read, q_start):
    """Return the query sequence from q_start to the end (includes soft-clipped bases)."""
    seq = read.query_sequence
    if seq is None or q_start is None:
        return ""
    if q_start < 0 or q_start >= len(seq):
        return ""
    return seq[q_start:]

def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences from alignments (primary + supplemental) starting at the reference start coordinate."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM/CRAM (indexed)")
    parser.add_argument("-c", "--chrom", required=True, help="Reference contig name")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start coordinate (1-based)")
    parser.add_argument("-l", "--length", type=int, required=True, help="Minimum length to keep")
    parser.add_argument("--sample", help="Optional sample name to prepend to read IDs")
    args = parser.parse_args()

    start_0 = args.start - 1
    required_len = args.length

    bam = pysam.AlignmentFile(args.bam, "rb")
    sample_name = args.sample if args.sample else get_sample_name(bam)

    seen_reads = set()  # track read names already output

    for read in bam.fetch(args.chrom, start_0, start_0 + 1):
        if read.is_unmapped or read.is_secondary:
            continue  # keep supplemental, skip secondary

        if read.query_name in seen_reads:
            continue  # skip duplicates

        qstart = find_qstart_for_ref_start(read, start_0)
        if qstart is None:
            continue

        seq = collect_sequence_from_qstart_to_end(read, qstart)
        if not seq or len(seq) < required_len:
            continue

        seen_reads.add(read.query_name)
        rid = f"{sample_name}_{read.query_name}" if sample_name else read.query_name
        print(f">{rid}_{len(seq)}\n{seq}")

    bam.close()

if __name__ == "__main__":
    main()
