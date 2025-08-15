#!/usr/bin/env python3
"""
aggregate_dotplot.py

Aggregate k-mer dot-matrices between a reference region FASTA and many "comp" FASTA files,
summing **presence across samples** (each (comp_pos, ref_pos) counted at most once per sample/file).

Produces PNG images of aggregated counts (total / forward /reverse) using scatter plots (one dot per
non-zero coordinate) and optionally saves numpy arrays (.npy).

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
    --palindrome-once \
    --marker-size 12
"""
import argparse
import glob
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
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

def parse_ref_coords_from_filename(reference_path):
    """
    Parse reference contig and start/end coordinate from the filename.
    Examples:
      chr16.161498.173978.fasta -> contig=chr16, start=161498, end=173978
      region_chr1_1000_2000.fa  -> take last two ints as start/end
    Returns (contig, start, end) or (None, None, None) if parsing fails.
    """
    stem = Path(reference_path).stem
    ints = re.findall(r"(\d+)", stem)
    if len(ints) >= 2:
        start = int(ints[-2])
        end = int(ints[-1])
        contig = stem.split(".")[0]
        return contig, start, end
    return None, None, None

def compute_max_comp_k(comp_fasta_paths, k):
    """
    Scan the comp FASTA files once to determine the maximum comp_k (comp_len - k + 1).
    Returns an integer >= 0.
    """
    paths = gather_paths(comp_fasta_paths)
    max_comp_k = 0
    for fasta_path in paths:
        for hdr, seq in iter_fasta_seqs(fasta_path):
            if seq is None:
                continue
            comp_k = len(seq) - k + 1
            if comp_k > max_comp_k:
                max_comp_k = comp_k
    return max_comp_k

# ---------- aggregation ----------
def aggregate_dot_matrices(reference_seq_file, comp_fasta_paths, k, palindrome_once=False, verbose=True):
    """
    Returns (forward_counts, reverse_counts, total_counts, ref_seq, comp_max_k)
    Matrices shape: (comp_max_k, ref_k) such that rows index comp k-mer starts (0-based),
    and columns index reference k-mer starts (0-based).

    Counting behavior:
      - For each comp FASTA file (treated as one "sample"), we mark whether a given
        (comp_pos, ref_pos) coordinate had any forward match, any reverse match, or both,
        across all sequences in that FASTA file.
      - After processing the entire FASTA file, we increment the aggregated matrices by 1
        for each coordinate that was observed in that sample in the appropriate matrix(es).
      - This ensures each sample contributes at most +1 to each matrix cell.
    """
    # --- Read reference ---
    with open(reference_seq_file, "r") as fh:
        ref_seq = "".join([l.strip() for l in fh if not l.startswith(">")]).upper()
    Lref = len(ref_seq)
    ref_k = Lref - k + 1
    if ref_k <= 0:
        raise ValueError(f"Reference region length {Lref} must be >= k ({k})")

    # --- Determine maximum comp k dimension across all comp FASTAs ---
    comp_max_k = compute_max_comp_k(comp_fasta_paths, k)
    if comp_max_k <= 0:
        raise ValueError("No comp sequences long enough for k across inputs")

    if verbose:
        print(f"Reference k-mers (ref_k) = {ref_k}; max comp k-mers (comp_max_k) = {comp_max_k}")

    # --- Initialize aggregated matrices with rows = comp_k, cols = ref_k ---
    forward_counts = np.zeros((comp_max_k, ref_k), dtype=np.int32)
    reverse_counts = np.zeros((comp_max_k, ref_k), dtype=np.int32)
    total_counts = np.zeros((comp_max_k, ref_k), dtype=np.int32)

    paths = gather_paths(comp_fasta_paths)
    if verbose:
        print(f"Found {len(paths)} FASTA files to process")

    path_iter = paths
    if HAVE_TQDM and verbose:
        path_iter = tqdm(paths, desc="FASTA files")

    for fasta_path in path_iter:
        # For this sample (fasta file) collect presence information:
        # dict mapping (j_comp, i_ref) -> {'f':bool, 'r':bool}
        sample_presence = {}

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

            mat = dp.mat  # mat orientation may be (ref_k, comp_k) or (comp_k, ref_k)

            # Convert sparse to CSR for efficient access if needed
            if hasattr(mat, "tocsc"):
                mat = mat.tocsr()

            # Normalize orientation so rows == ref_k and cols == comp_k.
            if mat.shape[0] != ref_k:
                if mat.shape[1] == ref_k:
                    mat = mat.T
                else:
                    # Attempt conservative trim
                    if mat.shape[0] > ref_k:
                        mat = mat[:ref_k, :]
                    elif mat.shape[1] > ref_k:
                        mat = mat[:, :ref_k]

            # Truncate columns to comp_max_k if necessary
            if mat.shape[1] > comp_max_k:
                mat = mat[:, :comp_max_k]

            # --- Find non-zero entries ---
            try:
                nz_rows, nz_cols = mat.nonzero()  # nz_rows: reference idx, nz_cols: comp idx
            except Exception:
                arr_mat = np.asarray(mat)
                nz = np.nonzero(arr_mat)
                nz_rows, nz_cols = nz[0], nz[1]

            # Force vals into flat numpy array of ints
            vals = np.array(mat[nz_rows, nz_cols]).ravel()

            # Update sample_presence for this FASTA file: mark forward/reverse presence per coordinate
            for i_ref, j_comp, v in zip(nz_rows, nz_cols, vals):
                # check bounds (defensive)
                if j_comp < 0 or j_comp >= comp_max_k or i_ref < 0 or i_ref >= ref_k:
                    continue
                coord = (j_comp, i_ref)
                state = sample_presence.get(coord)
                if state is None:
                    state = {"f": False, "r": False}
                    sample_presence[coord] = state

                if v == 2:  # palindrome
                    if palindrome_once:
                        state["f"] = True
                    else:
                        state["f"] = True
                        state["r"] = True
                elif v == 1:
                    state["f"] = True
                elif v == -1:
                    state["r"] = True

        # After processing all sequences in this FASTA file, update aggregated matrices
        for (j_comp, i_ref), state in sample_presence.items():
            any_seen = state["f"] or state["r"]
            if any_seen:
                total_counts[j_comp, i_ref] += 1
            if state["f"]:
                forward_counts[j_comp, i_ref] += 1
            if state["r"]:
                reverse_counts[j_comp, i_ref] += 1

    return forward_counts, reverse_counts, total_counts, ref_seq, comp_max_k

# ---------- plotting (scatter-only) ----------
def plot_aggregated_matrix(counts, output_png, x_coords=None, title=None, cmap='viridis',
                           vmax=None, logscale=False, label_prefix=None,
                           k=None, max_ticks=10, figsize=None, y_label_prefix="comp",
                           marker_size=20.0):
    """
    Plot scatter of counts with shape (comp_rows, ref_cols).
    - x_coords: iterable of reference start coordinates (length == ref_cols). If None, uses 1..ref_cols.
    - marker_size: matplotlib scatter 's' parameter (area in points^2)
    """
    arr = counts.astype(np.float64)
    if logscale:
        arr = np.log1p(arr)

    # select non-zero entries
    nz_mask = arr > 0
    if not nz_mask.any():
        # create an empty plot with titles/axes but no points
        nrows, ncols = arr.shape
        if figsize is None:
            width = min(max(4, ncols / 200 * 10), 20)
            height = min(max(4, nrows / 200 * 10), 20)
            figsize_use = (width, height)
        else:
            figsize_use = figsize
        fig, ax = plt.subplots(figsize=figsize_use)
        ax.set_title(label_prefix + " - " + (title or "Aggregated dot-matrix"))
        ax.set_xlabel('reference k-mer start (1-based)' if x_coords is None else 'reference base (k-mer start, 1-based)')
        ax.set_ylabel(f'{y_label_prefix} k-mer start (1-based)')
        plt.tight_layout()
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {output_png} (empty - no non-zero entries)")
        return

    ys, xs = np.nonzero(nz_mask)  # ys: comp idx (row), xs: ref idx (col)
    values = arr[ys, xs]

    # determine vmax robustly if not provided (use 99th percentile of values)
    if vmax is None:
        try:
            auto_vmax = float(np.percentile(values, 99))
            if auto_vmax <= 0:
                auto_vmax = float(values.max())
            if auto_vmax <= 0:
                auto_vmax = 1.0
        except Exception:
            auto_vmax = float(values.max()) if values.size else 1.0
        vmax_use = auto_vmax
    else:
        vmax_use = vmax

    cmap_obj = plt.get_cmap(cmap)
    # ensure colormap can handle "bad" values; not needed for scatter but keep robust behavior
    try:
        cmap_with_transparency = cmap_obj.copy()
    except Exception:
        try:
            colors = cmap_obj(np.linspace(0, 1, 256))
            cmap_with_transparency = mpl.colors.ListedColormap(colors)
        except Exception:
            cmap_with_transparency = cmap_obj

    if hasattr(cmap_with_transparency, "set_bad"):
        cmap_with_transparency.set_bad(alpha=0.0)

    nrows, ncols = arr.shape
    if figsize is None:
        width = min(max(4, ncols / 200 * 10), 20)
        height = min(max(4, nrows / 200 * 10), 20)
        figsize_use = (width, height)
    else:
        figsize_use = figsize

    fig, ax = plt.subplots(figsize=figsize_use)

    # scatter: x is reference index (col), y is comp index (row)
    sc = ax.scatter(xs, ys, c=values, cmap=cmap_with_transparency, s=marker_size, vmax=vmax_use, marker='s')
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label('log1p(counts)' if logscale else 'counts')

    if title is None:
        title = "Aggregated dot-matrix"
    if label_prefix:
        ax.set_title(f"{label_prefix} - {title}")
    else:
        ax.set_title(title)

    # tick placement helper
    def choose_ticks(n, max_ticks):
        if n <= max_ticks:
            return np.arange(n)
        locs = np.linspace(0, n - 1, num=max_ticks, dtype=int)
        return np.unique(locs)

    xticks = choose_ticks(ncols, max_ticks)
    yticks = choose_ticks(nrows, max_ticks)

    # X tick labels: use x_coords (reference base starts) if provided, else 1-based indices
    if x_coords is not None:
        x_coords_arr = np.asarray(x_coords)
        if x_coords_arr.shape[0] != ncols:
            xtick_labels = [str(int(x) + 1) for x in xticks]
        else:
            xtick_labels = [str(int(x_coords_arr[x])) for x in xticks]
        ax.set_xlabel('reference base (k-mer start, 1-based)')
    else:
        xtick_labels = [str(int(x) + 1) for x in xticks]
        ax.set_xlabel('reference k-mer start (1-based)')

    # Y tick labels: comp k-mer start positions (1-based)
    ytick_labels = [str(int(y) + 1) for y in yticks]
    ax.set_ylabel(f'{y_label_prefix} k-mer start (1-based)')

    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels, fontsize=8)

    # Invert y-axis to match imshow with origin='lower' semantics if desired.
    # The scatter used here has origin at (0,0) bottom-left due to plotting indices as-is.
    # Keep as-is so (0,0) is bottom-left similar to previous imshow(origin='lower') behavior.

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_png}")

# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser(description="Aggregate k-mer dotplots across many comp FASTA files (scatter plots).")
    parser.add_argument("--reference", required=True, help="Reference FASTA (one sequence covering start+length region)")
    parser.add_argument("--comps", required=True,
                        help="Comma-separated list or glob of comp FASTA files (each may be multi-FASTA). "
                             "Example: 'reads_dir/*.fasta' or 'a.fasta,b.fasta'")
    parser.add_argument("--k", type=int, required=True, help="kmer size")
    parser.add_argument("--output-root", required=True, help="Root prefix for output files (PNG and .npy)")
    parser.add_argument("--mode", default="total,forward,reverse",
                        help="Which matrices to produce/plot. Comma-separated subset of: total,forward,reverse")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
    parser.add_argument("--vmax", type=float, default=None, help="Optional vmax for scatter color scale")
    parser.add_argument("--logscale", action="store_true", help="Plot on log(1+x) scale")
    parser.add_argument("--threshold", type=int, default=None, help="Optional threshold to mask low counts (for plotting only)")
    parser.add_argument("--palindrome-once", action="store_true",
                        help="Count palindromic matches (value==2) only once into forward (instead of both forward+reverse)")
    parser.add_argument("--no-npy", action="store_true", help="Do not save .npy files of aggregated matrices")
    parser.add_argument("--progress", action="store_true", help="Show tqdm progress bar when available")
    parser.add_argument("--marker-size", type=float, default=0.01, help="Marker size (s) for scatter plot, in points^2")
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

    forward_counts, reverse_counts, total_counts, ref_seq, comp_max_k = aggregate_dot_matrices(
        args.reference, args.comps, args.k, palindrome_once=args.palindrome_once, verbose=True)

    # derive x_coords (reference base starts) from filename if possible
    contig, ref_start, ref_end = parse_ref_coords_from_filename(args.reference)
    ref_k = len(ref_seq) - args.k + 1
    if ref_start is not None:
        x_coords = np.arange(ref_start, ref_start + ref_k)
    else:
        x_coords = np.arange(1, ref_k + 1)  # 1-based positions fallback

    # Optionally apply threshold for plotting only (doesn't change saved arrays)
    def maybe_mask(arr):
        if args.threshold is not None:
            a = arr.copy()
            a[a < args.threshold] = 0
            return a
        return arr

    outroot = args.output_root
    if "total" in modes:
        arr = maybe_mask(total_counts)
        if not args.no_npy:
            np.save(f"{outroot}.total.npy", total_counts)
        plot_aggregated_matrix(arr, f"{outroot}.total.png",
                               x_coords=x_coords,
                               title=f"Total matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Total",
                               k=args.k, max_ticks=10, y_label_prefix="comp", marker_size=args.marker_size)

    if "forward" in modes:
        arr = maybe_mask(forward_counts)
        if not args.no_npy:
            np.save(f"{outroot}.forward.npy", forward_counts)
        plot_aggregated_matrix(arr, f"{outroot}.forward.png",
                               x_coords=x_coords,
                               title=f"Forward matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Forward",
                               k=args.k, max_ticks=10, y_label_prefix="comp", marker_size=args.marker_size)

    if "reverse" in modes:
        arr = maybe_mask(reverse_counts)
        if not args.no_npy:
            np.save(f"{outroot}.reverse.npy", reverse_counts)
        plot_aggregated_matrix(arr, f"{outroot}.reverse.png",
                               x_coords=x_coords,
                               title=f"Reverse matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Reverse",
                               k=args.k, max_ticks=10, y_label_prefix="comp", marker_size=args.marker_size)

    print("Done.")

if __name__ == "__main__":
    main()
