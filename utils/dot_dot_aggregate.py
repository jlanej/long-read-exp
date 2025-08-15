#!/usr/bin/env python3
"""
aggregate_dotplot.py

Aggregate k-mer dot-matrices between a reference region FASTA and many "comp" FASTA files,
summing **presence across samples** (each (comp_bin, ref_bin) counted at most once per sample/file).

Produces PNG images of aggregated counts (total / forward / reverse) using scatter plots
(one dot per non-zero binned coordinate) and optionally saves numpy arrays (.npy).

New CLI options:
  --alpha FLOAT          # marker alpha/transparency (0..1)
  --jitter FLOAT         # jitter fraction in bin units (0..1); spreads overlapping points
  --size-by-value        # scale marker size by the cell value
  --size-scale FLOAT     # multiplier for marker sizes when using --size-by-value

Binning:
  --bin-size N / --bin-ref-size / --bin-comp-size
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
    paths = []
    if isinstance(comps_arg, (list, tuple)):
        parts = list(comps_arg)
    elif isinstance(comps_arg, str):
        parts = [p.strip() for p in comps_arg.split(",")]
    else:
        parts = [comps_arg]

    for p in parts:
        if any(ch in p for ch in "*?[]"):
            paths.extend(sorted(glob.glob(p)))
        else:
            paths.append(p)

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
    stem = Path(reference_path).stem
    ints = re.findall(r"(\d+)", stem)
    if len(ints) >= 2:
        start = int(ints[-2])
        end = int(ints[-1])
        contig = stem.split(".")[0]
        return contig, start, end
    return None, None, None

def compute_max_comp_k(comp_fasta_paths, k):
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
    with open(reference_seq_file, "r") as fh:
        ref_seq = "".join([l.strip() for l in fh if not l.startswith(">")]).upper()
    Lref = len(ref_seq)
    ref_k = Lref - k + 1
    if ref_k <= 0:
        raise ValueError(f"Reference region length {Lref} must be >= k ({k})")

    comp_max_k = compute_max_comp_k(comp_fasta_paths, k)
    if comp_max_k <= 0:
        raise ValueError("No comp sequences long enough for k across inputs")

    if verbose:
        print(f"Reference k-mers (ref_k) = {ref_k}; max comp k-mers (comp_max_k) = {comp_max_k}")
        print(f"Binning: ref bin size = {bin_size_ref}; comp bin size = {bin_size_comp}")

    ref_bins = (ref_k + bin_size_ref - 1) // bin_size_ref
    comp_bins = (comp_max_k + bin_size_comp - 1) // bin_size_comp

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
        sample_presence_bins = {}

        for hdr, comp_seq in iter_fasta_seqs(fasta_path):
            if comp_seq is None:
                continue

            comp_len = len(comp_seq)
            comp_k = comp_len - k + 1
            if comp_k <= 0:
                continue

            try:
                dp = wotplot.DotPlotMatrix(ref_seq, comp_seq, k, binary=False, yorder="TB", verbose=False)
            except Exception as e:
                print(f"Error computing DotPlotMatrix for {fasta_path}:{hdr} -> {e}", file=sys.stderr)
                continue

            mat = dp.mat
            if hasattr(mat, "tocsc"):
                mat = mat.tocsr()

            if mat.shape[0] != ref_k:
                if mat.shape[1] == ref_k:
                    mat = mat.T
                else:
                    if mat.shape[0] > ref_k:
                        mat = mat[:ref_k, :]
                    elif mat.shape[1] > ref_k:
                        mat = mat[:, :ref_k]

            if mat.shape[1] > comp_max_k:
                mat = mat[:, :comp_max_k]

            try:
                nz_rows, nz_cols = mat.nonzero()
            except Exception:
                arr_mat = np.asarray(mat)
                nz = np.nonzero(arr_mat)
                nz_rows, nz_cols = nz[0], nz[1]

            vals = np.array(mat[nz_rows, nz_cols]).ravel()

            for i_ref, j_comp, v in zip(nz_rows, nz_cols, vals):
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
                           marker_size=20.0, alpha=0.9, jitter=0.0,
                           size_by_value=False, size_scale=1.0, rng_seed=None):
    """
    Scatter plot of binned counts.

    jitter: fraction of a bin (0..1) to jitter points by in both axes to reduce overplotting.
    size_by_value: if True, scale marker size by value (sqrt scaling).
    alpha: marker transparency (0..1).
    """
    arr = counts.astype(np.float64)
    if logscale:
        arr = np.log1p(arr)

    nz_mask = arr > 0
    if not nz_mask.any():
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

    ys, xs = np.nonzero(nz_mask)  # bin indices
    values = arr[ys, xs]

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

    # jitter (in bin units). small random offsets to reduce overlap
    rng = np.random.default_rng(rng_seed)
    if jitter and jitter > 0.0:
        jitter_x = (rng.random(xs.shape) - 0.5) * jitter
        jitter_y = (rng.random(ys.shape) - 0.5) * jitter
        xs_plot = xs + jitter_x
        ys_plot = ys + jitter_y
    else:
        xs_plot = xs.astype(float)
        ys_plot = ys.astype(float)

    # optional size scaling by value
    if size_by_value:
        vmax_val = values.max() if values.size else 1.0
        # sqrt scaling to compress dynamic range
        rel = (values / vmax_val)
        sizes = (np.sqrt(rel) * marker_size * size_scale) + (marker_size * 0.1)
    else:
        sizes = np.full_like(values, fill_value=marker_size, dtype=float)

    sc = ax.scatter(xs_plot, ys_plot, c=values, cmap=cmap_with_transparency,
                    s=sizes, vmax=vmax_use, alpha=alpha, marker='s', edgecolors='none')

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label('log1p(counts)' if logscale else 'counts')

    if title is None:
        title = "Aggregated dot-matrix"
    if label_prefix:
        ax.set_title(f"{label_prefix} - {title}")
    else:
        ax.set_title(title)

    def choose_ticks(n, max_ticks):
        if n <= max_ticks:
            return np.arange(n)
        locs = np.linspace(0, n - 1, num=max_ticks, dtype=int)
        return np.unique(locs)

    xticks = choose_ticks(ncols, max_ticks)
    yticks = choose_ticks(nrows, max_ticks)

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
    parser.add_argument("--reference", required=True)
    parser.add_argument("--comps", required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--mode", default="total,forward,reverse")
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--logscale", action="store_true")
    parser.add_argument("--threshold", type=int, default=None)
    parser.add_argument("--palindrome-once", action="store_true")
    parser.add_argument("--no-npy", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--marker-size", type=float, default=20.0, help="Base marker size in points^2")
    parser.add_argument("--alpha", type=float, default=0.9, help="Marker alpha (0..1)")
    parser.add_argument("--jitter", type=float, default=0.0, help="Jitter fraction of bin (0..1) to spread overlapping points")
    parser.add_argument("--size-by-value", action="store_true", help="Scale marker size by value (sqrt scaling)")
    parser.add_argument("--size-scale", type=float, default=1.0, help="Multiplier when using --size-by-value")
    parser.add_argument("--bin-size", type=int, default=None)
    parser.add_argument("--bin-ref-size", type=int, default=None)
    parser.add_argument("--bin-comp-size", type=int, default=None)
    args = parser.parse_args()

    global HAVE_TQDM
    if args.progress:
        HAVE_TQDM = True

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

    contig, ref_start, ref_end = parse_ref_coords_from_filename(args.reference)
    ref_k = len(ref_seq) - args.k + 1

    if ref_start is not None:
        raw_coords = np.arange(ref_start, ref_start + ref_k)
        bin_starts = np.arange(0, ref_k, bin_size_ref)
        bin_centers = []
        for s in bin_starts:
            end = min(s + bin_size_ref - 1, ref_k - 1)
            center_idx = s + (end - s) // 2
            bin_centers.append(raw_coords[center_idx])
        x_coords = np.array(bin_centers)
    else:
        x_coords = np.arange(1, ref_bins + 1)

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
                               k=args.k, max_ticks=10, y_label_prefix=f"comp (bin size={bin_size_comp})",
                               marker_size=args.marker_size, alpha=args.alpha, jitter=args.jitter,
                               size_by_value=args.size_by_value, size_scale=args.size_scale)

    if "forward" in modes:
        arr = maybe_mask(forward_counts)
        if not args.no_npy:
            np.save(f"{outroot}.forward.npy", forward_counts)
        plot_aggregated_matrix(arr, f"{outroot}.forward.png",
                               x_coords=x_coords,
                               title=f"Forward matches (k={args.k}, bin_ref={bin_size_ref}, bin_comp={bin_size_comp}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Forward",
                               k=args.k, max_ticks=10, y_label_prefix=f"comp (bin size={bin_size_comp})",
                               marker_size=args.marker_size, alpha=args.alpha, jitter=args.jitter,
                               size_by_value=args.size_by_value, size_scale=args.size_scale)

    if "reverse" in modes:
        arr = maybe_mask(reverse_counts)
        if not args.no_npy:
            np.save(f"{outroot}.reverse.npy", reverse_counts)
        plot_aggregated_matrix(arr, f"{outroot}.reverse.png",
                               x_coords=x_coords,
                               title=f"Reverse matches (k={args.k}, bin_ref={bin_size_ref}, bin_comp={bin_size_comp}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Reverse",
                               k=args.k, max_ticks=10, y_label_prefix=f"comp (bin size={bin_size_comp})",
                               marker_size=args.marker_size, alpha=args.alpha, jitter=args.jitter,
                               size_by_value=args.size_by_value, size_scale=args.size_scale)

    print("Done.")

if __name__ == "__main__":
    main()
