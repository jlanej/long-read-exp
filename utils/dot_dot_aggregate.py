#!/usr/bin/env python3
"""
aggregate_dotplot.py

Aggregate k-mer dot-matrices between a reference region FASTA and many "comp" FASTA files,
summing **presence across samples** (each (comp_bin, ref_bin) counted at most once per sample/file).

Produces PNG images of aggregated counts (total / forward / reverse) using scatter plots
(one dot per non-zero binned coordinate) and optionally saves numpy arrays (.npy).

New CLI options (binning):
  --bin-size N         # same bin size for both axes
  --bin-ref-size N     # bin size for reference axis only (overrides --bin-size for ref)
  --bin-comp-size N    # bin size for comp axis only (overrides --bin-size for comp)

Dependencies:
  - numpy
  - matplotlib
  - wotplot (https://github.com/fedarko/wotplot)
  - optionally tqdm (for progress bar)

Install with:
  pip install numpy matplotlib wotplot tqdm

Example:
  python aggregate_dotplot.py \
    --reference ref_region.fasta \
    --comps "reads_dir/*.fasta" \
    --k 11 \
    --output-root ref_region_agg \
    --mode total,forward,reverse \
    --logscale \
    --palindrome-once \
    --marker-size 12 \
    --bin-size 5
"""
import argparse
import glob
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
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

# ---------- aggregation with binning ----------
def aggregate_dot_matrices(reference_seq_file, comp_fasta_paths, k,
                           palindrome_once=False, verbose=True,
                           bin_size_ref=1, bin_size_comp=1):
    """
    Returns (forward_counts, reverse_counts, total_counts, ref_seq, comp_bins)
    Matrices shape: (comp_bins, ref_bins), where bins are non-overlapping windows
    of width bin_size_comp (rows) and bin_size_ref (cols) over k-mer start indices.

    Counting behavior:
      - Each FASTA file is treated as a sample. For each sample we collect presence
        of matches per raw coordinate, then map coordinates to bins, and then
        increment each binned cell by at most +1 per sample (per forward/reverse/total).
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
        print(f"Binning: ref bin size = {bin_size_ref}; comp bin size = {bin_size_comp}")

    # compute number of bins (ceil division)
    ref_bins = (ref_k + bin_size_ref - 1) // bin_size_ref
    comp_bins = (comp_max_k + bin_size_comp - 1) // bin_size_comp

    # --- Initialize aggregated binned matrices with rows = comp_bins, cols = ref_bins ---
    forward_counts = np.zeros((comp_bins, ref_bins), dtype=np.int32)
    reverse_counts = np.zeros((comp_bins, ref_bins), dtype=np.int32)
    total_counts = np.zeros((comp_bins, ref_bins), dtype=np.int32)

    paths = gather_paths(comp_fasta_paths)
    if verbose:
        print(f"Found {len(paths)} FASTA files to process")

    path_iter = paths
    if HAVE_TQDM and verbose:
        path_iter = tqdm(paths, desc="FASTA files")

    for fasta_path in path_iter:
        # per-sample binned presence map: (j_bin, i_bin) -> {'f':bool, 'r':bool}
        sample_presence_bins = {}

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

            vals = np.array(mat[nz_rows, nz_cols]).ravel()

            # Map raw coords to bins and update sample_presence_bins
            for i_ref, j_comp, v in zip(nz_rows, nz_cols, vals):
                # defensive bounds
                if j_comp < 0 or j_comp >= comp_max_k or i_ref < 0 or i_ref >= ref_k:
                    continue
                j_bin = j_comp // bin_size_comp
                i_bin = i_ref // bin_size_ref
                key = (j_bin, i_bin)
                state = sample_presence_bins.get(key)
                if state is None:
                    state = {"f": False, "r": False}
                    sample_presence_bins[key] = state

                if v == 2:
                    if palindrome_once:
                        state["f"] = True
                    else:
                        state["f"] = True
                        state["r"] = True
                elif v == 1:
                    state["f"] = True
                elif v == -1:
                    state["r"] = True

        # After processing the whole FASTA file, add +1 for each observed bin-state
        for (j_bin, i_bin), state in sample_presence_bins.items():
            if state["f"] or state["r"]:
                total_counts[j_bin, i_bin] += 1
            if state["f"]:
                forward_counts[j_bin, i_bin] += 1
            if state["r"]:
                reverse_counts[j_bin, i_bin] += 1

    return forward_counts, reverse_counts, total_counts, ref_seq, comp_bins, ref_bins

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
        # empty plot
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

    ys, xs = np.nonzero(nz_mask)  # ys: comp bin idx (row), xs: ref bin idx (col)
    values = arr[ys, xs]

    # determine vmax robustly if not provided (99th percentile)
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

    # X tick labels: use x_coords (reference bin centers) if provided, else 1-based bin indices
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

    ytick_labels = [str(int(y) + 1) for y in yticks]
    ax.set_ylabel(f'{y_label_prefix} k-mer start (1-based)')

    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_png}")

# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser(description="Aggregate k-mer dotplots across many comp FASTA files (scatter plots) with optional binning.")
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
    parser.add_argument("--marker-size", type=float, default=20.0, help="Marker size (s) for scatter plot, in points^2")
    parser.add_argument("--bin-size", type=int, default=None, help="Bin size applied to both ref and comp axes (overridden by axis-specific flags)")
    parser.add_argument("--bin-ref-size", type=int, default=None, help="Bin size for reference axis (k-mer starts)")
    parser.add_argument("--bin-comp-size", type=int, default=None, help="Bin size for comp axis (k-mer starts)")
    args = parser.parse_args()

    # respect --progress preference if tqdm is present
    global HAVE_TQDM
    if args.progress:
        HAVE_TQDM = True

    # determine bin sizes (default 1 => no binning)
    bin_size_ref = 1
    bin_size_comp = 1
    if args.bin_size is not None:
        if args.bin_size < 1:
            parser.error("--bin-size must be >= 1")
        bin_size_ref = args.bin_size
        bin_size_comp = args.bin_size
    if args.bin_ref_size is not None:
        if args.bin_ref_size < 1:
            parser.error("--bin-ref-size must be >= 1")
        bin_size_ref = args.bin_ref_size
    if args.bin_comp_size is not None:
        if args.bin_comp_size < 1:
            parser.error("--bin-comp-size must be >= 1")
        bin_size_comp = args.bin_comp_size

    modes = [m.strip().lower() for m in args.mode.split(",") if m.strip()]
    valid_modes = {"total", "forward", "reverse"}
    for m in modes:
        if m not in valid_modes:
            parser.error(f"Invalid mode: {m}. Choose subset of total,forward,reverse")

    forward_counts, reverse_counts, total_counts, ref_seq, comp_bins, ref_bins = aggregate_dot_matrices(
        args.reference, args.comps, args.k, palindrome_once=args.palindrome_once, verbose=True,
        bin_size_ref=bin_size_ref, bin_size_comp=bin_size_comp)

    # derive x_coords (reference bin centers) from filename if possible
    contig, ref_start, ref_end = parse_ref_coords_from_filename(args.reference)
    ref_k = len(ref_seq) - args.k + 1
    # compute bin centers (1-based coordinates)
    if ref_start is not None:
        # raw k-mer start coordinates (0-based): ref_start .. ref_start + ref_k - 1
        raw_coords = np.arange(ref_start, ref_start + ref_k)
        # take first element of each bin, then compute center
        bin_starts = np.arange(0, ref_k, bin_size_ref)
        bin_centers = []
        for s in bin_starts:
            # center in raw coordinates, using integer center
            end = min(s + bin_size_ref - 1, ref_k - 1)
            center_idx = s + (end - s) // 2
            bin_centers.append(raw_coords[center_idx])
        x_coords = np.array(bin_centers)
    else:
        # no absolute coords; supply 1-based bin indices
        x_coords = np.arange(1, ref_bins + 1)

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
                               title=f"Total matches (k={args.k}, bin_ref={bin_size_ref}, bin_comp={bin_size_comp}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Total",
                               k=args.k, max_ticks=10, y_label_prefix=f"comp (bin size={bin_size_comp})", marker_size=args.marker_size)

    if "forward" in modes:
        arr = maybe_mask(forward_counts)
        if not args.no_npy:
            np.save(f"{outroot}.forward.npy", forward_counts)
        plot_aggregated_matrix(arr, f"{outroot}.forward.png",
                               x_coords=x_coords,
                               title=f"Forward matches (k={args.k}, bin_ref={bin_size_ref}, bin_comp={bin_size_comp}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Forward",
                               k=args.k, max_ticks=10, y_label_prefix=f"comp (bin size={bin_size_comp})", marker_size=args.marker_size)

    if "reverse" in modes:
        arr = maybe_mask(reverse_counts)
        if not args.no_npy:
            np.save(f"{outroot}.reverse.npy", reverse_counts)
        plot_aggregated_matrix(arr, f"{outroot}.reverse.png",
                               x_coords=x_coords,
                               title=f"Reverse matches (k={args.k}, bin_ref={bin_size_ref}, bin_comp={bin_size_comp}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Reverse",
                               k=args.k, max_ticks=10, y_label_prefix=f"comp (bin size={bin_size_comp})", marker_size=args.marker_size)

    print("Done.")

if __name__ == "__main__":
    main()
