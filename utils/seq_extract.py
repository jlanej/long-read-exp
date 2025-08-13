#!/usr/bin/env python3
import pysam
import argparse

def get_sample_name(bam):
    """Extract sample name (SM) from BAM header if available."""
    header = bam.header
    samples = set()
    for rg in header.get("RG", []):
        if "SM" in rg:
            samples.add(rg["SM"])
    if len(samples) == 1:
        return samples.pop()
    elif len(samples) > 1:
        raise ValueError(f"Multiple sample names found in BAM header: {samples}")
    else:
        return None

def get_query_coords_for_region(read, start, end):
    """
    Given a read and a reference region [start, end) in 0-based coordinates,
    return the corresponding query sequence coordinates (q_start, q_end) in the read.
    Includes insertions inside the region.
    """
    ref_pos = read.reference_start
    query_pos = 0
    q_start = None
    q_end = None

    for (cigar_op, length) in read.cigartuples:
        if cigar_op == 0:  # match/mismatch
            for _ in range(length):
                if ref_pos == start:
                    q_start = query_pos
                if ref_pos == end:
                    q_end = query_pos
                    return q_start, q_end
                ref_pos += 1
                query_pos += 1
        elif cigar_op == 1:  # insertion (in query only)
            if q_start is not None and q_end is None:
                query_pos += length
            else:
                query_pos += length
        elif cigar_op == 2 or cigar_op == 3:  # deletion or skip (in reference only)
            for _ in range(length):
                if ref_pos == start:
                    q_start = query_pos
                if ref_pos == end:
                    q_end = query_pos
                    return q_start, q_end
                ref_pos += 1
        elif cigar_op in (4, 5):  # soft/hard clip
            query_pos += length

    # Handle case where region ends exactly at read end
    if q_start is not None and q_end is None:
        q_end = query_pos
    return q_start, q_end

def main():
    parser = argparse.ArgumentParser(
        description="Extract raw clipped sequences from BAM for a given reference range."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (indexed)")
    parser.add_argument("-c", "--chrom", required=True, help="Reference sequence/contig name")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start coordinate (1-based)")
    parser.add_argument("-e", "--end", type=int, required=True, help="End coordinate (1-based, inclusive)")
    parser.add_argument("--sample", help="Optional sample name to append to read IDs")
    args = parser.parse_args()

    # Convert start/end to 0-based, half-open interval
    start_0 = args.start - 1
    end_0 = args.end

    bam = pysam.AlignmentFile(args.bam, "rb")
    sample_name = args.sample if args.sample else get_sample_name(bam)

    for read in bam.fetch(args.chrom, start_0, end_0):
        # Require read to fully span the region
        if read.reference_start > start_0 or read.reference_end < end_0:
            continue

        q_start, q_end = get_query_coords_for_region(read, start_0, end_0)
        if q_start is None or q_end is None:
            continue  # couldn't map properly

        subseq = read.query_sequence[q_start:q_end]
        if subseq:
            if sample_name:
                read_id = f"{sample_name}_{read.query_name}"
            else:
                read_id = read.query_name
            print(f">{read_id}\n{subseq}")

    bam.close()

if __name__ == "__main__":
    main()
