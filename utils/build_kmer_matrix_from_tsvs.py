#!/usr/bin/env python3
"""
build_kmer_matrix_from_tsvs.py

Merge per-sample kmer TSVs (kmer<TAB>count) into a single kmer x sample matrix.

Usage:
  python build_kmer_matrix_from_tsvs.py samples/*.tsv --out kmers_matrix.tsv.gz --min-total-count 2 --min-samples-seen 10

Notes:
 - Each input TSV should be two columns: kmer <TAB> count (no header). Counts can be integer or float.
 - Filenames are used to derive sample names (basename without .tsv/.tsv.gz/.txt/.gz).
 - Output's first column will be 'kmer' and subsequent columns one sample each.
"""
from __future__ import annotations
import argparse
import glob
import os
from typing import List
import pandas as pd
from tqdm import tqdm

def expand_inputs(patterns: List[str]) -> List[str]:
    files = []
    for p in patterns:
        files.extend(sorted(glob.glob(p)))
    return files

def infer_sample_name(path: str) -> str:
    base = os.path.basename(path)
    for ext in (".tsv.gz", ".tsv", ".txt.gz", ".txt", ".kmer.gz", ".kmer", ".gz"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base

def load_single_tsv(path: str, sample_name: str) -> pd.Series:
    """
    Load one kmer TSV into a pandas Series (index=kmer, value=count).
    Accepts gzipped files; expects two whitespace-separated columns (kmer count).
    """
    s = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["kmer", sample_name],
        dtype={0: str, 1: "float64"},
        compression="infer",
        engine="python",
    )
    s = s.set_index("kmer")[sample_name]
    # drop explicit zero-count lines to reduce memory while merging; missing kmers will be interpreted as 0 later
    s = s[s != 0]
    return s

def build_matrix(file_list: List[str], min_total_count: int = 1, min_samples_seen: int = 1) -> pd.DataFrame:
    """
    Build matrix from a list of TSV paths. Returns DataFrame with rows=kmer, cols=samples.
    Applies filters:
      - keep kmers with total count across samples >= min_total_count
      - keep kmers seen (count>0) in at least min_samples_seen samples
    """
    series_list = []
    sample_order = []
    for path in tqdm(file_list, desc="Loading TSVs"):
        sample = infer_sample_name(path)
        sample_order.append(sample)
        s = load_single_tsv(path, sample)
        series_list.append(s)

    if not series_list:
        raise SystemExit("No input files loaded.")

    # outer-join all series into a dataframe; resulting df has rows=kmer, cols=samples
    df = pd.concat(series_list, axis=1).fillna(0)

    # ensure columns are in the order of input files
    df = df.loc[:, sample_order]

    # convert to integer if all values are integer-like
    try:
        # only attempt if no NaNs
        arr = df.values
        if (arr == arr.astype(int)).all():
            df = df.astype("int64")
    except Exception:
        # leave as float if conversion not safe
        pass

    # apply filters
    if min_total_count > 1:
        tot = df.sum(axis=1)
        df = df.loc[tot >= min_total_count]

    if min_samples_seen > 1:
        seen = (df > 0).sum(axis=1)
        df = df.loc[seen >= min_samples_seen]

    # move index name
    df.index.name = "kmer"
    return df

def main():
    p = argparse.ArgumentParser(description="Merge per-sample kmer TSVs into a kmer x sample matrix.")
    p.add_argument("inputs", nargs="+", help="Input TSV paths or glob patterns (e.g. results/*.tsv)")
    p.add_argument("--out", required=True, help="Output TSV path (can end with .gz for compression)")
    p.add_argument("--min-total-count", type=int, default=1,
                   help="Drop kmers with total count across samples < this (default: 1)")
    p.add_argument("--min-samples-seen", type=int, default=10,
                   help="Drop kmers seen (count>0) in fewer than this many samples (default: 10)")
    args = p.parse_args()

    files = expand_inputs(args.inputs)
    if not files:
        raise SystemExit("No input files matched the given patterns.")

    print(f"[INFO] {len(files)} files to process. Examples: {files[:5]}")
    df = build_matrix(files, min_total_count=args.min_total_count, min_samples_seen=args.min_samples_seen)
    print(f"[INFO] Final matrix: {df.shape[0]} kmers x {df.shape[1]} samples")

    # write out (allow .gz)
    compression = "gzip" if args.out.endswith(".gz") else None
    df.to_csv(args.out, sep="\t", compression=compression)
    print(f"[INFO] Wrote matrix to {args.out}")

if __name__ == "__main__":
    main()
