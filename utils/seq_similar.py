#!/usr/bin/env python3
"""
jaccard_upgma.py

Compute exact k-mer Jaccard distances (all-vs-all) from a multi-FASTA
and build a UPGMA tree (Newick).

Outputs (in outdir):
  - distances.tsv       : square distance matrix (tab-separated)
  - upgma_tree.newick   : Newick tree (UPGMA)

Requirements:
  - Python 3.8+
  - biopython
  - pandas
  - numpy
  - scipy

Usage example:
  python jaccard_upgma.py -i asm.all.fasta -o results_folder --k 21 --threads 4
"""
import argparse
import os
import sys
import time
from Bio import SeqIO
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from Bio.Phylo.Newick import Clade, Tree
from Bio import Phylo
from multiprocessing import Pool

# Global shared kmers list for multiprocessing workers (set by initializer)
GLOBAL_KMER_LIST = None

def build_kmer_set(seq, k):
    """Return a set of k-mers from seq (uppercase). Skips k-mers containing 'N'."""
    s = seq.upper()
    n = len(s)
    if n < k:
        return frozenset()
    km = set()
    end = n - k + 1
    for i in range(end):
        kmer = s[i:i+k]
        if 'N' in kmer:
            continue
        km.add(kmer)
    return frozenset(km)

def worker_pair(args):
    """Worker: given pair (i,j) compute Jaccard distance using GLOBAL_KMER_LIST."""
    i, j = args
    Ai = GLOBAL_KMER_LIST[i]
    Aj = GLOBAL_KMER_LIST[j]
    # both empty -> identical (distance 0)
    if not Ai and not Aj:
        return (i, j, 0.0)
    # iterate over smaller set for intersection
    if len(Ai) < len(Aj):
        inter = sum(1 for x in Ai if x in Aj)
    else:
        inter = sum(1 for x in Aj if x in Ai)
    union = len(Ai) + len(Aj) - inter
    if union == 0:
        jaccard = 1.0
    else:
        jaccard = inter / union
    dist = 1.0 - jaccard
    return (i, j, float(dist))

def linkage_to_biopython_tree(Z, labels):
    node = to_tree(Z, rd=False)
    def build_clade(n):
        if n.is_leaf():
            return Clade(branch_length=n.dist, name=labels[n.id])
        left = build_clade(n.get_left())
        right = build_clade(n.get_right())
        return Clade(branch_length=n.dist, clades=[left,right])
    return Tree(root=build_clade(node))

def parse_fasta_build_kmers(fasta, k, show_progress=True):
    recs = list(SeqIO.parse(fasta, "fasta"))
    if len(recs) == 0:
        raise SystemExit(f"No sequences found in {fasta}")
    ids = [r.id for r in recs]
    seqs = [str(r.seq) for r in recs]
    n = len(ids)
    kmers = []
    t0 = time.time()
    for i, s in enumerate(seqs):
        if show_progress and (i % 50 == 0):
            print(f"  Building k-mers: {i}/{n} sequences (elapsed {time.time()-t0:.1f}s)", file=sys.stderr)
        km = build_kmer_set(s, k)
        kmers.append(km)
    return ids, kmers

def main():
    p = argparse.ArgumentParser(description="Exact k-mer Jaccard all-vs-all + UPGMA")
    p.add_argument("-i", "--input", required=True, help="Input multi-FASTA")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--k", type=int, default=31, help="k-mer size (default 31)")
    p.add_argument("--threads", type=int, default=3, help="Number of worker processes (default 3)")
    p.add_argument("--no-symmetric-check", action="store_true", help="Skip symmetric check (slightly faster)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Reading FASTA and building k-mer sets (k={args.k})...", file=sys.stderr)
    ids, kmers_list = parse_fasta_build_kmers(args.input, args.k)

    n = len(ids)
    print(f"Built k-mers for {n} sequences. Starting pairwise Jaccard computations...", file=sys.stderr)

    # prepare pair list (i < j)
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]

    # use multiprocessing Pool with global kmers (to avoid pickling huge sets repeatedly)
    global GLOBAL_KMER_LIST
    GLOBAL_KMER_LIST = kmers_list

    results = []
    t0 = time.time()
    if args.threads > 1:
        print(f"Using {args.threads} processes.", file=sys.stderr)
        with Pool(processes=args.threads, initializer=_init_worker, initargs=(kmers_list,)) as pool:
            # pool.map returns results in order of pairs
            for i, res in enumerate(pool.imap_unordered(worker_pair, pairs, chunksize=256), 1):
                results.append(res)
                if i % 500 == 0:
                    print(f"  Completed {i}/{len(pairs)} pairs (elapsed {time.time()-t0:.1f}s)", file=sys.stderr)
    else:
        # single-threaded
        for idx, pair in enumerate(pairs, 1):
            res = worker_pair(pair)
            results.append(res)
            if idx % 500 == 0:
                print(f"  Completed {idx}/{len(pairs)} pairs (elapsed {time.time()-t0:.1f}s)", file=sys.stderr)

    # Build full matrix and fill diag
    mat = np.zeros((n, n), dtype=float)
    for (i, j, dist) in results:
        mat[i, j] = dist
        mat[j, i] = dist
    np.fill_diagonal(mat, 0.0)

    # Save distance matrix as TSV with IDs
    df = pd.DataFrame(mat, index=ids, columns=ids)
    dist_file = os.path.join(args.outdir, "distances.tsv")
    df.to_csv(dist_file, sep="\t", index=True)
    print(f"Wrote distance matrix to {dist_file}", file=sys.stderr)

    # Ensure symmetry
    if not args.no_symmetric_check:
        if not np.allclose(df.values, df.values.T):
            print("Warning: matrix not symmetric; symmetrizing.", file=sys.stderr)
            M = (df.values + df.values.T) / 2.0
            df = pd.DataFrame(M, index=ids, columns=ids)

    # Build condensed array and UPGMA (average linkage)
    condensed = squareform(df.values)
    print("Building UPGMA tree (average linkage)...", file=sys.stderr)
    Z = linkage(condensed, method="average")
    tree = linkage_to_biopython_tree(Z, ids)
    newick_file = os.path.join(args.outdir, "upgma_tree.newick")
    Phylo.write(tree, newick_file, "newick")
    print(f"Wrote UPGMA tree to {newick_file}", file=sys.stderr)
    print("Done.", file=sys.stderr)

def _init_worker(kmers):
    # initializer for Pool: set global kmers list in worker processes to avoid repeated pickling
    global GLOBAL_KMER_LIST
    GLOBAL_KMER_LIST = kmers

if __name__ == "__main__":
    main()
