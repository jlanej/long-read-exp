#!/usr/bin/env python3
import argparse
import pysam
from tqdm import tqdm
import sys

def revcomp(s: str) -> str:
    """Reverse complement a DNA string."""
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def parse_region(region: str):
    """
    Parse UCSC-style region strings like 'chr1:1000-2000'.
    Returns (contig, start, end) or (contig, None, None).
    """
    if ":" not in region:
        return region, None, None
    chrom, coords = region.split(":")
    if "-" in coords:
        start, end = coords.split("-")
        return chrom, int(start), int(end)
    else:
        return chrom, int(coords), None

def write_with_mate(read, in_f, out_f, written, kmer, kmer_rc, mate_failures):
    """Check read for k-mer, write it and its mate if present."""
    seq = read.query_sequence
    if not seq:
        return

    seq = seq.upper()
    if kmer in seq or kmer_rc in seq:
        if read.query_name not in written:
            out_f.write(read)
            written.add(read.query_name)
        try:
            mate = in_f.mate(read)
            if mate and mate.query_name not in written:
                out_f.write(mate)
                written.add(mate.query_name)
        except ValueError as e:
            # mate not retrievable
            mate_failures.append((read.query_name, str(e)))

def filter_cram_by_kmer_rare(input_cram, output_cram, ref_fasta, kmer, region=None):
    kmer = kmer.upper()
    kmer_rc = revcomp(kmer)

    region_args = {}
    if region:
        chrom, start, end = parse_region(region)
        region_args = dict(contig=chrom, start=start, end=end)

    written = set()
    mate_failures = []

    with pysam.AlignmentFile(input_cram, "rc", reference_filename=ref_fasta) as in_f, \
         pysam.AlignmentFile(output_cram, "wc", header=in_f.header, reference_filename=ref_fasta) as out_f:

        iterator = in_f.fetch(**region_args, until_eof=(not region))
        for read in tqdm(iterator, desc="Scanning CRAM"):
            write_with_mate(read, in_f, out_f, written, kmer, kmer_rc, mate_failures)

    # Reporting
    if mate_failures:
        sys.stderr.write(
            f"\n[WARNING] Could not retrieve mates for {len(mate_failures)} reads:\n"
        )
        for qname, reason in mate_failures[:20]:  # cap to avoid flooding
            sys.stderr.write(f"  {qname}: {reason}\n")
        if len(mate_failures) > 20:
            sys.stderr.write(f"  ... {len(mate_failures) - 20} more not shown\n")
    else:
        sys.stderr.write("\n[INFO] All mates retrieved successfully.\n")

def main():
    parser = argparse.ArgumentParser(
        description="Extract reads from CRAM containing a given k-mer (optimized for rare kmers)."
    )
    parser.add_argument("input_cram", help="Input CRAM file")
    parser.add_argument("output_cram", help="Output CRAM file")
    parser.add_argument("ref_fasta", help="Reference FASTA file")
    parser.add_argument("kmer", help="K-mer to search for")
    parser.add_argument(
        "--region",
        help="Restrict to UCSC-style region (e.g. chr1:1000-2000). "
             "If omitted, scans whole file.",
    )
    args = parser.parse_args()

    filter_cram_by_kmer_rare(
        args.input_cram,
        args.output_cram,
        args.ref_fasta,
        args.kmer,
        region=args.region,
    )

if __name__ == "__main__":
    main()
