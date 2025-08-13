#!/usr/bin/env python3
import pysam
import argparse
from collections import defaultdict

def get_sample_name(bam):
    samples = {rg.get("SM") for rg in bam.header.get("RG", []) if "SM" in rg}
    samples.discard(None)
    if len(samples) == 1:
        return samples.pop()
    elif len(samples) > 1:
        raise ValueError(f"Multiple sample names in BAM header: {samples}")
    return None

def merge_intervals(intervals):
    """Merge sorted intervals (start,end) and return merged list."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:  # overlap or touch
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged

def intervals_cover_region(merged_intervals, start_0, end_0):
    """Return True if union of merged_intervals covers [start_0, end_0) without gaps."""
    if not merged_intervals:
        return False
    # merged_intervals assumed sorted and non-overlapping
    # check starts at or before start_0 and ends at or after end_0 and no internal gap
    s0, e0 = merged_intervals[0]
    if s0 > start_0:
        return False
    cur_end = e0
    for s, e in merged_intervals[1:]:
        if s > cur_end:  # gap
            return False
        cur_end = max(cur_end, e)
    return cur_end >= end_0

def collect_aligned_pairs_for_read(alignments):
    """
    Collect aligned_pairs from all alignments for a read.
    Return dictionary qpos -> refpos (prefers non-None refpos if multiple entries).
    """
    q_to_r = {}
    for aln in alignments:
        try:
            pairs = aln.get_aligned_pairs(matches_only=False, with_seq=False)
        except ValueError:
            pairs = []
        for qpos, rpos in pairs:
            if qpos is None:
                # this is a pure deletion reference-only position - nothing to map to query
                continue
            if qpos not in q_to_r:
                q_to_r[qpos] = rpos
            else:
                # prefer non-None refpos over None
                if q_to_r[qpos] is None and rpos is not None:
                    q_to_r[qpos] = rpos
                # if both non-None and different, keep the first (rare)
    return q_to_r

def reconstruct_subsequence_from_qmap(q_to_r, query_sequence, start_0, end_0):
    """
    q_to_r: dict mapping query index -> ref index (ref index may be None)
    Return subsequence built in query order:
      - include bases where ref in [start, end)
      - include insertion bases (ref None) if they are located between two ref positions with at least one inside [start,end)
    """
    if not q_to_r:
        return ""

    # Build sorted list of (qpos, rpos) in ascending qpos
    items = sorted(q_to_r.items(), key=lambda x: x[0])
    qpos_list = [q for q, r in items]
    rpos_list = [r for q, r in items]

    # Precompute nearest previous non-none ref index for each position
    prev_ref = [None] * len(items)
    last = None
    for i, r in enumerate(rpos_list):
        if r is not None:
            last = r
        prev_ref[i] = last

    # Precompute nearest next non-none ref index for each position
    next_ref = [None] * len(items)
    nxt = None
    for i in range(len(items)-1, -1, -1):
        if rpos_list[i] is not None:
            nxt = rpos_list[i]
        next_ref[i] = nxt

    subseq_chars = []
    L = len(query_sequence)
    for i, (qpos, rpos) in enumerate(items):
        # safety check qpos in bounds
        if qpos is None or qpos < 0 or qpos >= L:
            continue
        include = False
        if rpos is not None:
            if start_0 <= rpos < end_0:
                include = True
        else:
            # insertion: include if either previous or next non-none ref pos falls within region
            p = prev_ref[i]
            n = next_ref[i]
            if (p is not None and start_0 <= p < end_0) or (n is not None and start_0 <= n < end_0):
                include = True
        if include:
            subseq_chars.append(query_sequence[qpos])
    return "".join(subseq_chars)

def extract_region_reads_merged(bam, chrom, start_0, end_0, sample_name):
    """Merge primary + supplementary alignments by read, ensure coverage, and reconstruct sequences."""
    # collect alignments that touch the region
    reads_dict = defaultdict(list)
    for read in bam.fetch(chrom, start_0, end_0):
        if read.is_unmapped or read.is_secondary:
            continue
        # include primary and supplementary; store alignment objects
        reads_dict[read.query_name].append(read)

    for rname, alignments in reads_dict.items():
        # compute merged reference intervals from the alignments
        ref_intervals = [(a.reference_start, a.reference_end) for a in alignments if a.reference_start is not None]
        if not ref_intervals:
            continue
        merged = merge_intervals(ref_intervals)
        if not intervals_cover_region(merged, start_0, end_0):
            # this read's alignments do not collectively span the region (there's a gap), skip
            continue

        # build qpos->refpos map across all alignments (query positions are global to read)
        q_to_r = collect_aligned_pairs_for_read(alignments)

        # reconstruct subsequence in query order, including insertions that sit inside region
        # need query_sequence from any alignment (all alignments refer to same read)
        query_sequence = None
        for a in alignments:
            if a.query_sequence:
                query_sequence = a.query_sequence
                break
        if not query_sequence:
            continue

        subseq = reconstruct_subsequence_from_qmap(q_to_r, query_sequence, start_0, end_0)
        if not subseq:
            continue

        rid = f"{sample_name}_{rname}" if sample_name else rname
        yield f"{rid}_{len(subseq)}", subseq

def main():
    parser = argparse.ArgumentParser(
        description="Extract merged read sequences spanning a reference region (merge supplementary alignments)."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM/CRAM (indexed)")
    parser.add_argument("-c", "--chrom", required=True, help="Reference contig name")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start coordinate (1-based)")
    parser.add_argument("-e", "--end", type=int, required=True, help="End coordinate (1-based, inclusive)")
    parser.add_argument("--sample", help="Optional sample name (overrides BAM header)")
    args = parser.parse_args()

    start_0 = args.start - 1
    end_0 = args.end

    bam = pysam.AlignmentFile(args.bam, "rb")
    sample_name = args.sample if args.sample else get_sample_name(bam)

    for read_id, seq in extract_region_reads_merged(bam, args.chrom, start_0, end_0, sample_name):
        print(f">{read_id}\n{seq}")

    bam.close()

if __name__ == "__main__":
    main()
