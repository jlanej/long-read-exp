#!/usr/bin/env python3
"""
fasta_clustering_pipeline.py

Pipeline to:
 1. Combine multiple FASTA files (and optional reference) into a multi-FASTA
 2. Run MSA (mafft / clustalo) if available
 3. Compute pairwise identity distance matrix from alignment
 4. Hierarchically cluster sequences and produce dendrogram + heatmap
 5. Detect large deletions relative to a reference sequence from the MSA and
    produce a simple plot of fraction-gapped columns (visualizes deletions)

Usage example:
  python fasta_clustering_pipeline.py -i ./fastas -r reference.fasta -o results --msa mafft

Dependencies (Python):
  pip install biopython scipy matplotlib numpy

External tools (optional, for MSA step if you want faster/better aligners):
  mafft or clustalo (the script will try mafft then clustalo)

Outputs (in the --outdir):
  - combined.fasta
  - aligned.fasta
  - distance_matrix.csv
  - dendrogram.png
  - heatmap_distance.png
  - gap_fraction_by_ref_position.png
  - deletion_intervals.bed

"""

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from collections import namedtuple

import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from scipy.cluster import hierarchy
from scipy.spatial import distance
import matplotlib.pyplot as plt


def check_program_exists(prog_name):
    return shutil.which(prog_name) is not None


def run_cmd(cmd, shell=False):
    print(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    res = subprocess.run(cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print("Command failed:")
        print(res.stderr)
        raise RuntimeError(f"Command failed: {cmd}")
    return res.stdout


def combine_fastas(input_dir, outfile, include_reference=None):
    paths = sorted(Path(input_dir).glob('*.fa')) + sorted(Path(input_dir).glob('*.fasta'))
    if include_reference:
        paths = [Path(include_reference)] + [p for p in paths if str(p) != str(include_reference)]
    if len(paths) == 0:
        raise FileNotFoundError(f"No fasta files found in {input_dir}")
    with open(outfile, 'w') as outfh:
        for p in paths:
            for rec in SeqIO.parse(p, 'fasta'):
                SeqIO.write(rec, outfh, 'fasta')
    print(f"Wrote combined FASTA with {len(list(SeqIO.parse(outfile,'fasta')))} records to {outfile}")


def run_msa(input_fasta, output_aln, prefer=None, threads=1):
    # if prefer provided, try that; otherwise try mafft then clustalo
    tried = []
    if prefer:
        tools = [prefer]
    else:
        tools = ['mafft', 'clustalo']
    for t in tools:
        if check_program_exists(t):
            print(f"Using MSA tool: {t}")
            if t == 'mafft':
                cmd = [t, '--auto', '--thread', str(threads), input_fasta]
                out = subprocess.run(cmd, stdout=open(output_aln, 'w'), stderr=subprocess.PIPE, text=True)
                if out.returncode == 0:
                    return output_aln
                else:
                    print(out.stderr)
            elif t == 'clustalo':
                cmd = [t, '-i', input_fasta, '-o', output_aln, '--force', '--outfmt=fasta']
                out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if out.returncode == 0:
                    return output_aln
                else:
                    print(out.stderr)
        else:
            tried.append(t)
    # If we didn't return, fail with helpful message
    raise RuntimeError(f"No supported MSA tool found or they failed. Tried: {tried}. Please install mafft or clustalo, or run your own MSA and pass aligned fasta to the script.")


def read_alignment(aln_path):
    # try fasta/aln formats
    return AlignIO.read(aln_path, 'fasta')


def compute_identity_distance_matrix(alignment: MultipleSeqAlignment):
    n = len(alignment)
    names = [rec.id for rec in alignment]
    seqs = [str(rec.seq) for rec in alignment]

    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            a = seqs[i]
            b = seqs[j]
            matches = 0
            valid_positions = 0
            for ca, cb in zip(a, b):
                # ignore columns where both are gaps
                if ca == '-' and cb == '-':
                    continue
                valid_positions += 1
                if ca == cb and ca != '-':
                    matches += 1
            if valid_positions == 0:
                identity = 0.0
            else:
                identity = matches / valid_positions
            dist = 1.0 - identity
            mat[i, j] = mat[j, i] = dist
    return names, mat


def write_distance_matrix(names, mat, outcsv):
    import csv
    with open(outcsv, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow([''] + names)
        for n, row in zip(names, mat):
            writer.writerow([n] + list(row))
    print(f"Saved distance matrix to {outcsv}")


def plot_dendrogram(mat, names, outpng, method='average'):
    # transform to condensed
    condensed = distance.squareform(mat)
    Z = hierarchy.linkage(condensed, method=method)
    plt.figure(figsize=(10, 6))
    dn = hierarchy.dendrogram(Z, labels=names, leaf_rotation=90)
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()
    print(f"Saved dendrogram to {outpng}")


def plot_heatmap(mat, names, outpng):
    plt.figure(figsize=(8, 8))
    plt.imshow(mat, aspect='auto', interpolation='nearest')
    plt.colorbar()
    plt.xticks(range(len(names)), names, rotation=90, fontsize=6)
    plt.yticks(range(len(names)), names, fontsize=6)
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()
    print(f"Saved heatmap to {outpng}")


def detect_deletions_vs_reference(alignment: MultipleSeqAlignment, reference_id=None, gap_fraction_threshold=0.5, min_length=100):
    # If reference_id given, try to find it. Otherwise assume first record is reference.
    names = [rec.id for rec in alignment]
    if reference_id and reference_id in names:
        ref_idx = names.index(reference_id)
    else:
        ref_idx = 0
        reference_id = names[0]

    ref_seq = str(alignment[ref_idx].seq)
    ncols = alignment.get_alignment_length()
    nseqs = len(alignment)

    gap_frac = np.zeros(ncols, dtype=float)
    for c in range(ncols):
        col = [str(alignment[i].seq)[c] for i in range(nseqs)]
        gaps = sum(1 for x in col if x == '-')
        gap_frac[c] = gaps / nseqs

    # Map alignment column indices to reference genome coordinates (1-based)
    ref_pos_by_col = []
    ref_pos = 0
    for c in range(ncols):
        if ref_seq[c] != '-':
            ref_pos += 1
            ref_pos_by_col.append(ref_pos)
        else:
            ref_pos_by_col.append(None)

    # Find contiguous regions where gap_frac > threshold and reference is non-gap
    intervals = []  # (ref_start, ref_end, mean_gap_fraction, aln_start_col, aln_end_col)
    in_interval = False
    start_col = None
    for c in range(ncols):
        if gap_frac[c] >= gap_fraction_threshold and ref_pos_by_col[c] is not None:
            if not in_interval:
                in_interval = True
                start_col = c
        else:
            if in_interval:
                end_col = c - 1
                # convert to ref coords
                ref_start = ref_pos_by_col[start_col]
                # find last non-None in the block
                ref_end = ref_pos_by_col[end_col]
                if ref_start is not None and ref_end is not None and (ref_end - ref_start + 1) >= min_length:
                    mean_gap = float(np.mean(gap_frac[start_col:end_col+1]))
                    intervals.append((ref_start, ref_end, mean_gap, start_col, end_col))
                in_interval = False
                start_col = None
    # if still in interval at end
    if in_interval and start_col is not None:
        end_col = ncols - 1
        ref_start = ref_pos_by_col[start_col]
        ref_end = ref_pos_by_col[end_col]
        if ref_start is not None and ref_end is not None and (ref_end - ref_start + 1) >= min_length:
            mean_gap = float(np.mean(gap_frac[start_col:end_col+1]))
            intervals.append((ref_start, ref_end, mean_gap, start_col, end_col))

    return gap_frac, ref_pos_by_col, intervals, reference_id


def plot_gap_fraction_vs_ref(gap_frac, ref_pos_by_col, outpng, reference_name):
    # Build arrays of reference positions (skip None)
    ref_positions = [p for p in ref_pos_by_col if p is not None]
    gap_values = [gap_frac[i] for i, p in enumerate(ref_pos_by_col) if p is not None]
    plt.figure(figsize=(10, 4))
    plt.plot(ref_positions, gap_values)
    plt.xlabel('Reference position')
    plt.ylabel('Fraction of sequences with gap at column')
    plt.title(f'Gap fraction vs reference position (ref={reference_name})')
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()
    print(f"Saved gap-fraction plot to {outpng}")


def save_intervals_bed(intervals, outbed, ref_name='ref'):
    with open(outbed, 'w') as fh:
        for (s, e, mean_gap, a, b) in intervals:
            # BED is 0-based half-open; convert 1-based ref coords to 0-based
            fh.write(f"{ref_name}\t{s-1}\t{e}\tdeletion;mean_gap={mean_gap:.3f};aln_cols={a}-{b}\n")
    print(f"Saved intervals to {outbed}")


def cluster_assignments_from_distance(mat, names, threshold=0.1):
    # Use linkage and fcluster
    condensed = distance.squareform(mat)
    Z = hierarchy.linkage(condensed, method='average')
    labels = hierarchy.fcluster(Z, t=threshold, criterion='distance')
    # return mapping name -> label
    return dict(zip(names, labels))


def main():
    parser = argparse.ArgumentParser(description='Cluster FASTA sequences and visualize deletions vs reference using MSA')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('-r', '--reference', required=False, help='Path to reference fasta (optional). If given, included first in combined FASTA')
    parser.add_argument('-o', '--outdir', default='results', help='Output directory')
    parser.add_argument('--msa', choices=['mafft', 'clustalo'], default=None, help='Preferred MSA tool; if not provided script tries available tools')
    parser.add_argument('--threads', type=int, default=1, help='Threads for MSA (if supported)')
    parser.add_argument('--gap-threshold', type=float, default=0.5, help='Fraction of sequences gapped at MSA column to call deletion region')
    parser.add_argument('--min-deletion-len', type=int, default=100, help='Minimum length in reference coordinates to report interval')
    parser.add_argument('--cluster-threshold', type=float, default=0.1, help='Distance threshold for cutting dendrogram into clusters')

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    combined = outdir / 'combined.fasta'
    aligned = outdir / 'aligned.fasta'
    dist_csv = outdir / 'distance_matrix.csv'
    dendro_png = outdir / 'dendrogram.png'
    heatmap_png = outdir / 'heatmap_distance.png'
    gap_plot = outdir / 'gap_fraction_by_ref_position.png'
    bed_out = outdir / 'deletion_intervals.bed'

    combine_fastas(args.input_dir, combined, include_reference=args.reference)

    # Run MSA (or expect user to have provided aligned fasta at outdir/aligned.fasta)
    try:
        run_msa(str(combined), str(aligned), prefer=args.msa, threads=args.threads)
    except Exception as e:
        print("MSA step failed or no external MSA tool available. If you already have an alignment, place it at results/aligned.fasta or run MSA separately and re-run this script.")
        raise

    alignment = read_alignment(str(aligned))

    names, mat = compute_identity_distance_matrix(alignment)
    write_distance_matrix(names, mat, dist_csv)

    plot_dendrogram(mat, names, dendro_png)
    plot_heatmap(mat, names, heatmap_png)

    # cluster assignments
    assignments = cluster_assignments_from_distance(mat, names, threshold=args.cluster_threshold)
    with open(outdir / 'cluster_assignments.tsv', 'w') as fh:
        fh.write('name\tcluster\n')
        for n in names:
            fh.write(f"{n}\t{assignments[n]}\n")
    print(f"Saved cluster assignments to {outdir / 'cluster_assignments.tsv'}")

    # Deletion detection vs reference
    gap_frac, ref_pos_by_col, intervals, ref_name = detect_deletions_vs_reference(alignment, reference_id=args.reference and Path(args.reference).stem or None, gap_fraction_threshold=args.gap_threshold, min_length=args.min_deletion_len)
    plot_gap_fraction_vs_ref(gap_frac, ref_pos_by_col, gap_plot, ref_name)
    save_intervals_bed(intervals, bed_out, ref_name)

    print('\nPipeline complete. Check the output directory for results.')


if __name__ == '__main__':
    main()
