#!/usr/bin/env python3
import argparse
import pysam
from collections import Counter
from tqdm import tqdm

def parse_region(region: str):
    """Parse UCSC-style region like 'chr1:1000-2000'."""
    if ":" not in region:
        return region, None, None
    chrom, coords = region.split(":")
    if "-" in coords:
        start, end = coords.split("-")
        return chrom, int(start), int(end)
    else:
        return chrom, int(coords), None

def revcomp(seq: str) -> str:
    """Reverse complement DNA sequence."""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def canonical_kmer(kmer: str) -> str:
    """Return canonical form of kmer = lexicographically smaller of kmer and revcomp."""
    kmer_rc = revcomp(kmer)
    return min(kmer, kmer_rc)

def count_kmers_in_cram(input_cram, ref_fasta, ksize, canonical=False, region=None):
    counts = Counter()

    region_args = {}
    if region:
        chrom, start, end = parse_region(region)
        region_args = dict(contig=chrom, start=start, end=end)

    with pysam.AlignmentFile(input_cram, "rc", reference_filename=ref_fasta) as in_f:
        iterator = in_f.fetch(**region_args, until_eof=(not region))
        for read in tqdm(iterator, desc="Scanning CRAM"):
            seq = read.query_sequence
            if not seq:
                continue
            seq = seq.upper()
            for i in range(len(seq) - ksize + 1):
                kmer = seq[i : i + ksize]
                if "N" in kmer:  # skip ambiguous kmers
                    continue
                if canonical:
                    kmer = canonical_kmer(kmer)
                counts[kmer] += 1

    return counts

def main():
    parser = argparse.ArgumentParser(
        description="Count k-mers in a CRAM file (optionally canonicalized, optionally restricted to a region)."
    )
    parser.add_argument("input_cram", help="Input CRAM file")
    parser.add_argument("ref_fasta", help="Reference FASTA file")
    parser.add_argument("ksize", type=int, help="K-mer size")
    parser.add_argument(
        "--canonical",
        action="store_true",
        help="Collapse k-mer and its reverse complement into one canonical count.",
    )
    parser.add_argument(
        "--region",
        help="Restrict to UCSC-style region (e.g. chr1:1000-2000). "
             "If omitted, scans the whole file.",
    )
    parser.add_argument(
        "--out",
        help="Output TSV file (kmer\\tcount). If omitted, prints to stdout.",
    )
    args = parser.parse_args()

    counts = count_kmers_in_cram(
        args.input_cram,
        args.ref_fasta,
        args.ksize,
        canonical=args.canonical,
        region=args.region,
    )

    # Output results
    if args.out:
        with open(args.out, "w") as fh:
            for kmer, count in counts.most_common():
                fh.write(f"{kmer}\t{count}\n")
    else:
        for kmer, count in counts.most_common():
            print(f"{kmer}\t{count}")

if __name__ == "__main__":
    main()
