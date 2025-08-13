#!/usr/bin/env python3
"""
aggregate_dotplot.py

Aggregate k-mer dot-matrices between a reference region FASTA and many "comp" FASTA files,
summing match counts across all comp sequences.

Produces PNG images of aggregated counts (total / forward / reverse) and optionally saves
numpy arrays (.npy).

Dependencies:
  - numpy
  - matplotlib
  - wotplot (https://github.com/fedarko/wotplot)
  - optionally tqdm (for progress bar)

Install with:
  pip install numpy matplotlib wotplot tqdm

Usage example:
  python aggregate_dotplot.py \
    --reference ref_region.fasta \
    --comps "reads_dir/*.fasta" \
    --k 11 \
    --output-root ref_region_agg \
    --mode total,forward,reverse \
    --logscale \
    --palindrome-once
"""

import argparse
import glob
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import wotplot
import sys

try:
    from tqdm import tqdm
    HAVE_TQDM = True
except Exception:
    HAVE_TQDM = False

# ---------- helpers ----------
def iter_fasta_seqs(fasta_path):
    """Yield (header, seq) for each sequence in a FASTA file (uppercased)."""
    hdr = None
    seq_chunks = []
    with open(fasta_path, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq_chunks).upper()
                hdr = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if hdr is not None:
            yield hdr, "".join(seq_chunks).upper()

def root_file_name_sans_dir(file_name):
    return Path(file_name).stem

# ---------- aggregation ----------
def gather_paths(comps_arg):
    """
    Accepts:
      - comma-separated list of filenames/globs
      - a single glob pattern
      - a list (if caller passes list)
    Returns ordered list of existing file paths.
    """
    paths = []
    if isinstance(comps_arg, (list, tuple)):
        parts = list(comps_arg)
    elif isinstance(comps_arg, str):
        # split on commas, allow globbing
        parts = [p.strip() for p in comps_arg.split(",")]
    else:
        parts = [comps_arg]

    for p in parts:
        if any(ch in p for ch in "*?[]"):
            paths.extend(sorted(glob.glob(p)))
        else:
            paths.append(p)

    # remove duplicates and only keep existing files
    seen = set()
    out = []
    for p in paths:
        if p in seen:
            continue
        seen.add(p)
        if Path(p).exists():
            out.append(p)
        else:
            print(f"Warning: path not found -> {p}", file=sys.stderr)
    return out

def aggregate_dot_matrices(reference_seq_file, comp_fasta_paths, k, palindrome_once=False, verbose=True):
    """
    Returns (forward_counts, reverse_counts, total_counts, ref_seq)
    """
    # --- Read reference ---
    with open(reference_seq_file, "r") as fh:
        ref_seq = "".join([l.strip() for l in fh if not l.startswith(">")]).upper()
    Lref = len(ref_seq)
    ref_k = Lref - k + 1
    if ref_k <= 0:
        raise ValueError(f"Reference region length {Lref} must be >= k ({k})")

    # --- Initialize aggregated matrices ---
    forward_counts = np.zeros((ref_k, ref_k), dtype=np.int32)
    reverse_counts = np.zeros((ref_k, ref_k), dtype=np.int32)
    total_counts = np.zeros((ref_k, ref_k), dtype=np.int32)

    paths = gather_paths(comp_fasta_paths)
    if verbose:
        print(f"Found {len(paths)} FASTA files to process")

    path_iter = paths
    if HAVE_TQDM and verbose:
        path_iter = tqdm(paths, desc="FASTA files")

    for fasta_path in path_iter:
        for hdr, comp_seq in iter_fasta_seqs(fasta_path):
            if comp_seq is None:
                continue

            comp_len = len(comp_seq)
            comp_k = comp_len - k + 1
            if comp_k <= 0:
                continue  # sequence too short for k-mers

            # --- Compute dotplot matrix ---
            try:
                dp = wotplot.DotPlotMatrix(ref_seq, comp_seq, k, binary=False, yorder="TB", verbose=False)
            except Exception as e:
                print(f"Error computing DotPlotMatrix for {fasta_path}:{hdr} -> {e}", file=sys.stderr)
                continue

            mat = dp.mat  # can be sparse or dense

            # Convert sparse to CSR for efficient row access
            if hasattr(mat, "tocsc"):
                mat = mat.tocsr()

            # Limit to reference length only (do not truncate comp side)
            mat = mat[:, :ref_k]

            # --- Find non-zero entries ---
            nz_rows, nz_cols = mat.nonzero()
            # Force vals into flat numpy array of ints
            vals = np.array(mat[nz_rows, nz_cols]).ravel()

            # --- Increment counts ---
            for i, j, v in zip(nz_rows, nz_cols, vals):
                total_counts[i, j] += 1
                if v == 2:  # palindrome
                    if palindrome_once:
                        forward_counts[i, j] += 1
                    else:
                        forward_counts[i, j] += 1
                        reverse_counts[i, j] += 1
                elif v == 1:
                    forward_counts[i, j] += 1
                elif v == -1:
                    reverse_counts[i, j] += 1

    return forward_counts, reverse_counts, total_counts, ref_seq


# ---------- plotting ----------
def plot_aggregated_matrix(counts, output_png, title=None, cmap='viridis', vmax=None, logscale=False, label_prefix=None):
    """
    Plot heatmap of counts (2D numpy array). Saves to output_png.
    """
    # choose display array
    arr = counts.astype(np.float64)
    if logscale:
        arr = np.log1p(arr)

    fig, ax = plt.subplots(figsize=(10,10))
    im = ax.imshow(arr, origin='lower', interpolation='nearest', cmap=cmap, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('log1p(counts)' if logscale else 'counts')
    if title is None:
        title = f"Aggregated dot-matrix"
    if label_prefix:
        title = f"{label_prefix} - {title}"
    ax.set_title(title)
    ax.set_xlabel('reference k-mer index (j)')
    ax.set_ylabel('reference k-mer index (i)')
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_png}")

# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser(description="Aggregate k-mer dotplots across many comp FASTA files.")
    parser.add_argument("--reference", required=True, help="Reference FASTA (one sequence covering start+length region)")
    parser.add_argument("--comps", required=True,
                        help="Comma-separated list or glob of comp FASTA files (each may be multi-FASTA). "
                             "Example: 'reads_dir/*.fasta' or 'a.fasta,b.fasta'")
    parser.add_argument("--k", type=int, required=True, help="kmer size")
    parser.add_argument("--output-root", required=True, help="Root prefix for output files (PNG and .npy)")
    parser.add_argument("--mode", default="total,forward,reverse",
                        help="Which matrices to produce/plot. Comma-separated subset of: total,forward,reverse")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
    parser.add_argument("--vmax", type=float, default=None, help="Optional vmax for imshow color scale")
    parser.add_argument("--logscale", action="store_true", help="Plot on log(1+x) scale")
    parser.add_argument("--threshold", type=int, default=None, help="Optional threshold to mask low counts (for plotting only)")
    parser.add_argument("--palindrome-once", action="store_true",
                        help="Count palindromic matches (value==2) only once into forward (instead of both forward+reverse)")
    parser.add_argument("--no-npy", action="store_true", help="Do not save .npy files of aggregated matrices")
    parser.add_argument("--progress", action="store_true", help="Show tqdm progress bar when available")
    args = parser.parse_args()

    # respect --progress preference if tqdm is present
    global HAVE_TQDM
    if args.progress:
        HAVE_TQDM = True

    modes = [m.strip().lower() for m in args.mode.split(",") if m.strip()]
    valid_modes = {"total", "forward", "reverse"}
    for m in modes:
        if m not in valid_modes:
            parser.error(f"Invalid mode: {m}. Choose subset of total,forward,reverse")

    forward_counts, reverse_counts, total_counts, ref_seq = aggregate_dot_matrices(
        args.reference, args.comps, args.k, palindrome_once=args.palindrome_once, verbose=True)

    # Optionally apply threshold for plotting only (doesn't change saved arrays)
    def maybe_mask(arr):
        if args.threshold is not None:
            mask = arr < args.threshold
            a = arr.copy()
            a[mask] = 0
            return a
        return arr

    # Save and plot requested modes
    outroot = args.output_root
    if "total" in modes:
        arr = maybe_mask(total_counts)
        if not args.no_npy:
            np.save(f"{outroot}.total.npy", total_counts)
        plot_aggregated_matrix(arr, f"{outroot}.total.png",
                               title=f"Total matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Total")

    if "forward" in modes:
        arr = maybe_mask(forward_counts)
        if not args.no_npy:
            np.save(f"{outroot}.forward.npy", forward_counts)
        plot_aggregated_matrix(arr, f"{outroot}.forward.png",
                               title=f"Forward matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Forward")

    if "reverse" in modes:
        arr = maybe_mask(reverse_counts)
        if not args.no_npy:
            np.save(f"{outroot}.reverse.npy", reverse_counts)
        plot_aggregated_matrix(arr, f"{outroot}.reverse.png",
                               title=f"Reverse matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Reverse")

    print("Done.")

if __name__ == "__main__":
    main()
