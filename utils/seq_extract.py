#!/usr/bin/env python3
import pysam
import argparse

def get_sample_name(bam):
    """Extract sample name (SM) from BAM header if available."""
    samples = {rg.get("SM") for rg in bam.header.get("RG", []) if "SM" in rg}
    samples.discard(None)
    if len(samples) == 1:
        return samples.pop()
    elif len(samples) > 1:
        raise ValueError(f"Multiple sample names in BAM header: {samples}")
    return None

def extract_region_reads(bam, chrom, start_0, end_0, sample_name, clip=False):
    for read in bam.fetch(chrom, start_0, end_0):
        if read.is_unmapped:
            continue

        # Must span region (or optionally clip to nearest aligned base)
        if not clip and (read.reference_start > start_0 or read.reference_end < end_0):
            continue

        aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)

        q_start = None
        q_end = None

        for qpos, rpos in aligned_pairs:
            # Find first query position >= start
            if q_start is None and rpos is not None and rpos >= start_0:
                q_start = qpos
            # Find first query position >= end
            if rpos is not None and rpos >= end_0:
                q_end = qpos
                break

        # Clip to nearest aligned base if start/end is in a gap
        if clip:
            if q_start is None:
                # Find first query pos after start in aligned_pairs
                for qpos, rpos in aligned_pairs:
                    if qpos is not None and (rpos is None or rpos > start_0):
                        q_start = qpos
                        break
            if q_end is None:
                for qpos, rpos in reversed(aligned_pairs):
                    if qpos is not None and (rpos is None or rpos < end_0):
                        q_end = qpos + 1  # include this base
                        break

        if q_start is None or q_end is None or q_end <= q_start:
            continue

        subseq = read.query_sequence[q_start:q_end]
        rid = f"{sample_name}_{read.query_name}" if sample_name else read.query_name
        yield rid, subseq

def main():
    parser = argparse.ArgumentParser(
        description="Extract raw clipped read sequences spanning a reference region."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM/CRAM (indexed)")
    parser.add_argument("-c", "--chrom", required=True, help="Reference contig name")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start coordinate (1-based)")
    parser.add_argument("-e", "--end", type=int, required=True, help="End coordinate (1-based, inclusive)")
    parser.add_argument("--sample", help="Optional sample name (overrides BAM header)")
    parser.add_argument("--clip", action="store_true",
                        help="Clip to nearest aligned base if start or end falls in a deletion/gap")
    args = parser.parse_args()

    start_0 = args.start - 1
    end_0 = args.end

    bam = pysam.AlignmentFile(args.bam, "rb")
    sample_name = args.sample if args.sample else get_sample_name(bam)

    for read_id, seq in extract_region_reads(bam, args.chrom, start_0, end_0, sample_name, clip=args.clip):
        print(f">{read_id}\n{seq}")

    bam.close()

if __name__ == "__main__":
    main()
