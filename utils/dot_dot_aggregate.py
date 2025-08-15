#!/usr/bin/env python3
"""
aggregate_dotplot.py

Aggregate k-mer dot-matrices between a reference region FASTA and many "comp" FASTA files,
summing presence across samples (each (comp_bin, ref_bin) counted at most once per sample/file).

Features:
 - scatter plotting (binned)
 - binning (--bin-size / --bin-ref-size / --bin-comp-size)
 - marker alpha, jitter, size-by-value
 - option to combine strands into a single signed matrix (--reverse-negative)
 - option to produce two panels with separate color scales (--split-panels)

Dependencies:
  - numpy
  - matplotlib
  - wotplot
  - optionally tqdm

Example:
  python aggregate_dotplot.py --reference ref.fasta --comps "reads/*.fasta" --k 11 \
      --output-root out --mode total --bin-size 5 --marker-size 16 --alpha 0.6 \
      --jitter 0.3 --size-by-value --reverse-negative
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
    """
    Returns (forward_counts, reverse_counts, total_counts, ref_seq, comp_bins, ref_bins)
    where each matrix is shape (comp_bins, ref_bins) containing integer sample counts.
    """
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
        # per-sample presence per bin
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

# ---------- plotting helpers (robust/sanitizing versions) ----------
from matplotlib.colors import TwoSlopeNorm

def scatter_plot_single(ax, xs, ys, values, *, cmap, vmax, logscale, marker_size,
                        alpha, jitter, size_by_value, size_scale, rng_seed=None):
    """
    Robust single-panel scatter helper that sanitizes input arrays before plotting.
    Returns the scatter object (or None if nothing plotted).
    """
    # convert to numpy arrays and make floats
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    values = np.asarray(values, dtype=float)

    # apply logscale to values if requested (sign-preserving)
    if logscale:
        # sign * log1p(abs)
        values = np.sign(values) * np.log1p(np.abs(values))

    # jitter in bin units if requested (apply only if shapes match)
    if jitter and jitter > 0.0 and xs.size == ys.size:
        rng = np.random.default_rng(rng_seed)
        jitter_x = (rng.random(xs.shape) - 0.5) * jitter
        jitter_y = (rng.random(ys.shape) - 0.5) * jitter
        xs_plot = xs + jitter_x
        ys_plot = ys + jitter_y
    else:
        xs_plot = xs.copy()
        ys_plot = ys.copy()

    # sanitize: drop any entries with non-finite coords or values
    finite_mask = np.isfinite(xs_plot) & np.isfinite(ys_plot) & np.isfinite(values)
    if not np.any(finite_mask):
        return None

    xs_plot = xs_plot[finite_mask]
    ys_plot = ys_plot[finite_mask]
    values = values[finite_mask]

    # compute sizes robustly
    if size_by_value:
        vmax_val = np.max(np.abs(values)) if values.size else 1.0
        if not np.isfinite(vmax_val) or vmax_val <= 0:
            vmax_val = 1.0
        rel = (np.abs(values) / vmax_val)
        sizes = (np.sqrt(rel) * marker_size * size_scale) + (marker_size * 0.1)
    else:
        sizes = np.full(values.shape, fill_value=marker_size, dtype=float)

    # sanitize sizes
    sizes = np.nan_to_num(sizes, nan=marker_size * 0.1, posinf=marker_size * 10.0, neginf=marker_size * 0.1)
    sizes[sizes <= 0] = marker_size * 0.1

    # determine vmax for color mapping robustly
    if vmax is None or (not np.isfinite(vmax)) or (vmax <= 0):
        try:
            vmax_use = float(np.percentile(np.abs(values), 99))
            if not np.isfinite(vmax_use) or vmax_use <= 0:
                vmax_use = float(np.max(np.abs(values))) if values.size else 1.0
        except Exception:
            vmax_use = float(np.max(np.abs(values))) if values.size else 1.0
        if not np.isfinite(vmax_use) or vmax_use <= 0:
            vmax_use = 1.0
    else:
        vmax_use = float(vmax)

    # scatter
    sc = ax.scatter(xs_plot, ys_plot, c=values, cmap=cmap,
                    s=sizes, vmax=vmax_use, alpha=alpha, marker='s', edgecolors='none')
    return sc


def plot_signed_matrix(counts_signed, output_png, x_coords=None, title=None,
                       cmap='RdBu_r', vmax=None, logscale=False, label_prefix=None,
                       max_ticks=10, figsize=None, y_label_prefix="comp",
                       marker_size=20.0, alpha=0.9, jitter=0.0,
                       size_by_value=False, size_scale=1.0, rng_seed=None):
    """
    Robust signed scatter plot: negative -> reverse, positive -> forward.
    Uses TwoSlopeNorm centered at 0. Sanitizes inputs and avoids NaN/Inf.
    """
    arr = np.asarray(counts_signed, dtype=float)

    # sign-preserving logscale
    if logscale:
        arr = np.sign(arr) * np.log1p(np.abs(arr))

    nz_mask = arr != 0
    nrows, ncols = arr.shape

    if not nz_mask.any():
        # produce an empty but labeled figure
        if figsize is None:
            width = min(max(4, ncols / 200 * 10), 20)
            height = min(max(4, nrows / 200 * 10), 20)
            figsize_use = (width, height)
        else:
            figsize_use = figsize
        fig, ax = plt.subplots(figsize=figsize_use)
        ax.set_title(label_prefix + " - " + (title or "Aggregated dot-matrix (signed)"))
        ax.set_xlabel('reference k-mer start (1-based)' if x_coords is None else 'reference base (k-mer start, 1-based)')
        ax.set_ylabel(f'{y_label_prefix} k-mer start (1-based)')
        plt.tight_layout()
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {output_png} (empty - no non-zero entries)")
        return

    ys, xs = np.nonzero(nz_mask)
    values = arr[ys, xs].astype(float)

    # sanitize finite entries
    finite_mask = np.isfinite(xs) & np.isfinite(ys) & np.isfinite(values)
    if not np.any(finite_mask):
        if figsize is None:
            width = min(max(4, ncols / 200 * 10), 20)
            height = min(max(4, nrows / 200 * 10), 20)
            figsize_use = (width, height)
        else:
            figsize_use = figsize
        fig, ax = plt.subplots(figsize=figsize_use)
        ax.set_title(label_prefix + " - " + (title or "Aggregated dot-matrix (signed)"))
        ax.set_xlabel('reference k-mer start (1-based)' if x_coords is None else 'reference base (k-mer start, 1-based)')
        ax.set_ylabel(f'{y_label_prefix} k-mer start (1-based)')
        plt.tight_layout()
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {output_png} (no finite entries)")
        return

    xs = xs[finite_mask].astype(float)
    ys = ys[finite_mask].astype(float)
    values = values[finite_mask]

    # determine symmetric vmax safely
    if vmax is None or (not np.isfinite(vmax)) or (vmax <= 0):
        vmax_use = np.max(np.abs(values)) if values.size else 1.0
        if not np.isfinite(vmax_use) or vmax_use <= 0:
            vmax_use = 1.0
    else:
        vmax_use = float(vmax)

    # compute sizes robustly
    if size_by_value:
        vmax_val = np.max(np.abs(values)) if values.size else 1.0
        if not np.isfinite(vmax_val) or vmax_val <= 0:
            vmax_val = 1.0
        rel = (np.abs(values) / vmax_val)
        sizes = (np.sqrt(rel) * marker_size * size_scale) + (marker_size * 0.1)
    else:
        sizes = np.full(values.shape, fill_value=marker_size, dtype=float)

    sizes = np.nan_to_num(sizes, nan=marker_size * 0.1, posinf=marker_size * 10.0, neginf=marker_size * 0.1)
    sizes[sizes <= 0] = marker_size * 0.1

    # prepare cmap and norm
    cmap_obj = plt.get_cmap(cmap)
    try:
        cmap_with_transparency = cmap_obj.copy()
    except Exception:
        colors = cmap_obj(np.linspace(0, 1, 256))
        cmap_with_transparency = mpl.colors.ListedColormap(colors)
    if hasattr(cmap_with_transparency, "set_bad"):
        cmap_with_transparency.set_bad(alpha=0.0)

    norm = TwoSlopeNorm(vmin=-vmax_use, vcenter=0.0, vmax=vmax_use)

    # plot
    fig, ax = plt.subplots(figsize=(min(max(4, ncols / 200 * 10), 20), min(max(4, nrows / 200 * 10), 20)))
    sc = ax.scatter(xs, ys, c=values, cmap=cmap_with_transparency, s=sizes, norm=norm,
                    alpha=alpha, marker='s', edgecolors='none')

    # colorbar via ScalarMappable so ticks are symmetric around zero
    sm = mpl.cm.ScalarMappable(cmap=cmap_with_transparency, norm=norm)
    sm.set_array(np.array([-vmax_use, 0.0, vmax_use]))
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('signed counts (forward positive; reverse negative)')

    title_use = title or "Aggregated dot-matrix (signed)"
    if label_prefix:
        ax.set_title(f"{label_prefix} - {title_use}")
    else:
        ax.set_title(title_use)

    def choose_ticks(n, max_ticks):
        if n <= max_ticks:
            return np.arange(n)
        locs = np.linspace(0, n - 1, num=max_ticks, dtype=int)
        return np.unique(locs)

    xticks = choose_ticks(ncols, max_ticks)
    yticks = choose_ticks(nrows, max_ticks)

    if x_coords is not None:
        x_coords_arr = np.asarray(x_coords)
        if x_coords_arr.shape[0] == ncols:
            xtick_labels = [str(int(x_coords_arr[x])) for x in xticks]
            ax.set_xlabel('reference base (k-mer start, 1-based)')
        else:
            xtick_labels = [str(int(x) + 1) for x in xticks]
            ax.set_xlabel('reference k-mer start (1-based)')
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


def plot_split_panels(forward_arr, reverse_arr, output_png, x_coords=None, title=None,
                      cmap_f='viridis', cmap_r='viridis', vmax_f=None, vmax_r=None,
                      logscale=False, label_prefix=None, max_ticks=10, figsize=None,
                      y_label_prefix="comp", marker_size=20.0, alpha=0.9, jitter=0.0,
                      size_by_value=False, size_scale=1.0, rng_seed=None):
    """
    Plot two side-by-side panels: left = forward counts, right = reverse counts.
    Each panel has its own colorbar (separate color scales).
    Robust to empty panels.
    """
    f_arr = np.asarray(forward_arr, dtype=float)
    r_arr = np.asarray(reverse_arr, dtype=float)
    if logscale:
        f_arr = np.log1p(f_arr)
        r_arr = np.log1p(r_arr)

    nrows, ncols = f_arr.shape
    if figsize is None:
        width = min(max(6, ncols / 200 * 14), 28)
        height = min(max(4, nrows / 200 * 10), 20)
        figsize_use = (width, height)
    else:
        figsize_use = figsize

    fig, (axl, axr) = plt.subplots(1, 2, figsize=figsize_use, gridspec_kw={'width_ratios': [1, 1]})

    # forward panel
    nzf = f_arr > 0
    if nzf.any():
        ys_f, xs_f = np.nonzero(nzf)
        vals_f = f_arr[ys_f, xs_f]
        if vmax_f is None:
            try:
                auto_vmax_f = float(np.percentile(vals_f, 99))
                if auto_vmax_f <= 0:
                    auto_vmax_f = float(vals_f.max())
                if auto_vmax_f <= 0:
                    auto_vmax_f = 1.0
            except Exception:
                auto_vmax_f = float(vals_f.max()) if vals_f.size else 1.0
            vmax_f_use = auto_vmax_f
        else:
            vmax_f_use = float(vmax_f)
        sc_f = scatter_plot_single(axl, xs_f, ys_f, vals_f, cmap=cmap_f, vmax=vmax_f_use,
                                   logscale=False, marker_size=marker_size, alpha=alpha, jitter=jitter,
                                   size_by_value=size_by_value, size_scale=size_scale, rng_seed=rng_seed)
        if sc_f is not None:
            cbar_f = fig.colorbar(sc_f, ax=axl)
            cbar_f.set_label('counts (forward)')
        else:
            axl.text(0.5, 0.5, 'no forward hits', ha='center', va='center')
    else:
        axl.text(0.5, 0.5, 'no forward hits', ha='center', va='center')

    # reverse panel (plot as positive but labeled reverse)
    nzr = r_arr > 0
    if nzr.any():
        ys_r, xs_r = np.nonzero(nzr)
        vals_r = r_arr[ys_r, xs_r]
        if vmax_r is None:
            try:
                auto_vmax_r = float(np.percentile(vals_r, 99))
                if auto_vmax_r <= 0:
                    auto_vmax_r = float(vals_r.max())
                if auto_vmax_r <= 0:
                    auto_vmax_r = 1.0
            except Exception:
                auto_vmax_r = float(vals_r.max()) if vals_r.size else 1.0
            vmax_r_use = auto_vmax_r
        else:
            vmax_r_use = float(vmax_r)
        sc_r = scatter_plot_single(axr, xs_r, ys_r, vals_r, cmap=cmap_r, vmax=vmax_r_use,
                                   logscale=False, marker_size=marker_size, alpha=alpha, jitter=jitter,
                                   size_by_value=size_by_value, size_scale=size_scale, rng_seed=rng_seed)
        if sc_r is not None:
            cbar_r = fig.colorbar(sc_r, ax=axr)
            cbar_r.set_label('counts (reverse)')
        else:
            axr.text(0.5, 0.5, 'no reverse hits', ha='center', va='center')
    else:
        axr.text(0.5, 0.5, 'no reverse hits', ha='center', va='center')

    title_use = title or "Forward (left) / Reverse (right)"
    if label_prefix:
        fig.suptitle(f"{label_prefix} - {title_use}")
    else:
        fig.suptitle(title_use)

    def choose_ticks(n, max_ticks):
        if n <= max_ticks:
            return np.arange(n)
        locs = np.linspace(0, n - 1, num=max_ticks, dtype=int)
        return np.unique(locs)

    xticks = choose_ticks(ncols, max_ticks)
    yticks = choose_ticks(nrows, max_ticks)

    # labels / ticks for both axes
    if x_coords is not None:
        x_coords_arr = np.asarray(x_coords)
        if x_coords_arr.shape[0] == ncols:
            xtick_labels = [str(int(x_coords_arr[x])) for x in xticks]
            xlabel = 'reference base (k-mer start, 1-based)'
        else:
            xtick_labels = [str(int(x) + 1) for x in xticks]
            xlabel = 'reference k-mer start (1-based)'
    else:
        xtick_labels = [str(int(x) + 1) for x in xticks]
        xlabel = 'reference k-mer start (1-based)'

    ytick_labels = [str(int(y) + 1) for y in yticks]

    for ax in (axl, axr):
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
        ax.set_yticks(yticks)
        ax.set_yticklabels(ytick_labels, fontsize=8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'{y_label_prefix} k-mer start (1-based)')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_png}")

# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser(description="Aggregate k-mer dotplots across many comp FASTA files (scatter plots) with options for strand visualization.")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--comps", required=True, help="Comma-separated list or glob of comp FASTA files")
    parser.add_argument("--k", type=int, required=True, help="kmer size")
    parser.add_argument("--output-root", required=True, help="Root prefix for output files (PNG and .npy)")
    parser.add_argument("--mode", default="total,forward,reverse", help="Which matrices to produce/plot. Comma-separated subset of: total,forward,reverse")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap (sequential) for single-strand panels")
    parser.add_argument("--diverging-cmap", default="RdBu_r", help="Diverging colormap for signed combined plot")
    parser.add_argument("--vmax", type=float, default=None, help="Optional vmax for scatter color scale (or abs vmax for signed)")
    parser.add_argument("--logscale", action="store_true", help="Plot on log(1+x) scale")
    parser.add_argument("--threshold", type=int, default=None, help="Optional threshold to mask low counts (for plotting only)")
    parser.add_argument("--palindrome-once", action="store_true", help="Count palindromic matches only once into forward")
    parser.add_argument("--no-npy", action="store_true", help="Do not save .npy files of aggregated matrices")
    parser.add_argument("--progress", action="store_true", help="Show tqdm progress bar when available")
    parser.add_argument("--marker-size", type=float, default=20.0, help="Base marker size in points^2")
    parser.add_argument("--alpha", type=float, default=0.9, help="Marker alpha (0..1)")
    parser.add_argument("--jitter", type=float, default=0.0, help="Jitter fraction of bin (0..1)")
    parser.add_argument("--size-by-value", action="store_true", help="Scale marker size by value (sqrt scaling)")
    parser.add_argument("--size-scale", type=float, default=1.0, help="Multiplier when using --size-by-value")
    parser.add_argument("--bin-size", type=int, default=None, help="Bin size applied to both axes")
    parser.add_argument("--bin-ref-size", type=int, default=None)
    parser.add_argument("--bin-comp-size", type=int, default=None)
    parser.add_argument("--reverse-negative", action="store_true", help="Plot reverse counts as negative values in the same signed plot (diverging colormap).")
    parser.add_argument("--split-panels", action="store_true", help="Produce two side-by-side panels with separate colorbars (forward vs reverse).")
    args = parser.parse_args()

    global HAVE_TQDM
    if args.progress:
        HAVE_TQDM = True

    # determine bin sizes
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

    # compute x_coords for bin centers when possible
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

    # If reverse-negative is requested, we will compute signed = forward - reverse
    # and plot that for whichever mode the user requested (commonly "total").
    if args.reverse_negative:
        # signed counts: forward - reverse
        signed = forward_counts.astype(np.int32) - reverse_counts.astype(np.int32)
        # Save signed .npy if desired
        if not args.no_npy:
            np.save(f"{outroot}.signed.npy", signed)
        # If user asked for split-panels as well, ignore split-panels (mutually exclusive)
        if args.split_panels:
            print("Note: --split-panels ignored when --reverse-negative is set.", file=sys.stderr)

        if "total" in modes:
            arr = maybe_mask(signed)
            vmax_use = args.vmax if args.vmax is not None else None
            plot_signed_matrix(arr, f"{outroot}.signed.png",
                               x_coords=x_coords,
                               title=f"Signed matches (forward positive, reverse negative; k={args.k}) - {root_file_name_sans_dir(args.reference)}",
                               cmap=args.diverging_cmap, vmax=vmax_use, logscale=args.logscale,
                               label_prefix="Signed", max_ticks=10, y_label_prefix="comp",
                               marker_size=args.marker_size, alpha=args.alpha, jitter=args.jitter,
                               size_by_value=args.size_by_value, size_scale=args.size_scale, rng_seed=None)

    else:
        # Not using signed combined plot. Respect modes (total/forward/reverse).
        if "total" in modes:
            arr = maybe_mask(total_counts)
            if not args.no_npy:
                np.save(f"{outroot}.total.npy", total_counts)
            plot_aggregated_matrix_params = dict(
                output_png=f"{outroot}.total.png",
                x_coords=x_coords,
                title=f"Total matches (k={args.k}, bin_ref={bin_size_ref}, bin_comp={bin_size_comp}) - {root_file_name_sans_dir(args.reference)}",
                cmap=args.cmap, vmax=args.vmax, logscale=args.logscale, label_prefix="Total",
                max_ticks=10, y_label_prefix=f"comp (bin size={bin_size_comp})",
                marker_size=args.marker_size, alpha=args.alpha, jitter=args.jitter,
                size_by_value=args.size_by_value, size_scale=args.size_scale
            )
            if args.split_panels:
                plot_split_panels(forward_counts, reverse_counts, **plot_aggregated_matrix_params)
            else:
                arr_plot = arr
                nrows, ncols = arr_plot.shape
                nz_mask = arr_plot > 0
                if not nz_mask.any():
                    fig, ax = plt.subplots(figsize=(6,4))
                    ax.set_title(plot_aggregated_matrix_params['title'])
                    ax.set_xlabel('reference k-mer start (1-based)' if x_coords is None else 'reference base (k-mer start, 1-based)')
                    ax.set_ylabel(plot_aggregated_matrix_params['y_label_prefix'])
                    plt.tight_layout()
                    plt.savefig(plot_aggregated_matrix_params['output_png'], dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print(f"Saved {plot_aggregated_matrix_params['output_png']} (empty - no non-zero entries)")
                else:
                    ys, xs = np.nonzero(nz_mask)
                    vals = arr_plot[ys, xs].astype(np.float64)
                    if args.logscale:
                        vals = np.log1p(vals)
                    if args.vmax is None:
                        try:
                            auto_vmax = float(np.percentile(vals, 99))
                            if auto_vmax <= 0:
                                auto_vmax = float(vals.max())
                            if auto_vmax <= 0:
                                auto_vmax = 1.0
                        except Exception:
                            auto_vmax = float(vals.max()) if vals.size else 1.0
                        vmax_use = auto_vmax
                    else:
                        vmax_use = args.vmax

                    cmap_obj = plt.get_cmap(args.cmap)
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

                    if args.size_by_value:
                        vmax_val = vals.max() if vals.size else 1.0
                        rel = (vals / vmax_val)
                        sizes = (np.sqrt(rel) * args.marker_size * args.size_scale) + (args.marker_size * 0.1)
                    else:
                        sizes = np.full_like(vals, fill_value=args.marker_size, dtype=float)

                    fig, ax = plt.subplots(figsize=(min(max(4, ncols/200*10),20), min(max(4, nrows/200*10),20)))
                    sc = ax.scatter(xs, ys, c=vals, cmap=cmap_with_transparency, s=sizes, vmax=vmax_use,
                                    alpha=args.alpha, marker='s', edgecolors='none')
                    cbar = fig.colorbar(sc, ax=ax)
                    cbar.set_label('log1p(counts)' if args.logscale else 'counts')
                    ax.set_title(plot_aggregated_matrix_params['title'])
                    def choose_ticks(n, max_ticks):
                        if n <= max_ticks:
                            return np.arange(n)
                        locs = np.linspace(0, n - 1, num=max_ticks, dtype=int)
                        return np.unique(locs)
                    xticks = choose_ticks(ncols, 10)
                    yticks = choose_ticks(nrows, 10)
                    if x_coords is not None:
                        x_coords_arr = np.asarray(x_coords)
                        if x_coords_arr.shape[0] == ncols:
                            xtick_labels = [str(int(x_coords_arr[x])) for x in xticks]
                            ax.set_xlabel('reference base (k-mer start, 1-based)')
                        else:
                            xtick_labels = [str(int(x) + 1) for x in xticks]
                            ax.set_xlabel('reference k-mer start (1-based)')
                    else:
                        xtick_labels = [str(int(x) + 1) for x in xticks]
                        ax.set_xlabel('reference k-mer start (1-based)')
                    ytick_labels = [str(int(y) + 1) for y in yticks]
                    ax.set_xticks(xticks)
                    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
                    ax.set_yticks(yticks)
                    ax.set_yticklabels(ytick_labels, fontsize=8)
                    ax.set_ylabel(f"comp (bin size={bin_size_comp}) k-mer start (1-based)")
                    plt.tight_layout()
                    plt.savefig(plot_aggregated_matrix_params['output_png'], dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print(f"Saved {plot_aggregated_matrix_params['output_png']}")

        if "forward" in modes:
            arr = maybe_mask(forward_counts)
            if not args.no_npy:
                np.save(f"{outroot}.forward.npy", forward_counts)
            nrows, ncols = arr.shape
            nz_mask = arr > 0
            if not nz_mask.any():
                fig, ax = plt.subplots(figsize=(6,4))
                ax.set_title(f"Forward matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}")
                ax.set_xlabel('reference k-mer start (1-based)' if x_coords is None else 'reference base (k-mer start, 1-based)')
                ax.set_ylabel(f'comp (bin size={bin_size_comp}) k-mer start (1-based)')
                plt.tight_layout()
                plt.savefig(f"{outroot}.forward.png", dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"Saved {outroot}.forward.png (empty - no non-zero entries)")
            else:
                ys_f, xs_f = np.nonzero(nz_mask)
                vals_f = arr[ys_f, xs_f].astype(np.float64)
                if args.logscale:
                    vals_f = np.log1p(vals_f)
                if args.vmax is None:
                    try:
                        auto_vmax = float(np.percentile(vals_f, 99))
                        if auto_vmax <= 0:
                            auto_vmax = float(vals_f.max())
                        if auto_vmax <= 0:
                            auto_vmax = 1.0
                    except Exception:
                        auto_vmax = float(vals_f.max()) if vals_f.size else 1.0
                    vmax_use = auto_vmax
                else:
                    vmax_use = args.vmax
                cmap_obj = plt.get_cmap(args.cmap)
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
                if args.size_by_value:
                    vmax_val = vals_f.max() if vals_f.size else 1.0
                    rel = (vals_f / vmax_val)
                    sizes = (np.sqrt(rel) * args.marker_size * args.size_scale) + (args.marker_size * 0.1)
                else:
                    sizes = np.full_like(vals_f, fill_value=args.marker_size, dtype=float)
                fig, ax = plt.subplots(figsize=(min(max(4, ncols/200*10),20), min(max(4, nrows/200*10),20)))
                sc = ax.scatter(xs_f, ys_f, c=vals_f, cmap=cmap_with_transparency, s=sizes, vmax=vmax_use,
                                alpha=args.alpha, marker='s', edgecolors='none')
                cbar = fig.colorbar(sc, ax=ax)
                cbar.set_label('log1p(counts)' if args.logscale else 'counts')
                ax.set_title(f"Forward matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}")
                def choose_ticks(n, max_ticks):
                    if n <= max_ticks:
                        return np.arange(n)
                    locs = np.linspace(0, n - 1, num=max_ticks, dtype=int)
                    return np.unique(locs)
                xticks = choose_ticks(ncols, 10)
                yticks = choose_ticks(nrows, 10)
                if x_coords is not None:
                    x_coords_arr = np.asarray(x_coords)
                    if x_coords_arr.shape[0] == ncols:
                        xtick_labels = [str(int(x_coords_arr[x])) for x in xticks]
                        ax.set_xlabel('reference base (k-mer start, 1-based)')
                    else:
                        xtick_labels = [str(int(x) + 1) for x in xticks]
                        ax.set_xlabel('reference k-mer start (1-based)')
                else:
                    xtick_labels = [str(int(x) + 1) for x in xticks]
                    ax.set_xlabel('reference k-mer start (1-based)')
                ytick_labels = [str(int(y) + 1) for y in yticks]
                ax.set_xticks(xticks)
                ax.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
                ax.set_yticks(yticks)
                ax.set_yticklabels(ytick_labels, fontsize=8)
                ax.set_ylabel(f"comp (bin size={bin_size_comp}) k-mer start (1-based)")
                plt.tight_layout()
                plt.savefig(f"{outroot}.forward.png", dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"Saved {outroot}.forward.png")

        if "reverse" in modes:
            arr = maybe_mask(reverse_counts)
            if not args.no_npy:
                np.save(f"{outroot}.reverse.npy", reverse_counts)
            nrows, ncols = arr.shape
            nz_mask = arr > 0
            if not nz_mask.any():
                fig, ax = plt.subplots(figsize=(6,4))
                ax.set_title(f"Reverse matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}")
                ax.set_xlabel('reference k-mer start (1-based)' if x_coords is None else 'reference base (k-mer start, 1-based)')
                ax.set_ylabel(f'comp (bin size={bin_size_comp}) k-mer start (1-based)')
                plt.tight_layout()
                plt.savefig(f"{outroot}.reverse.png", dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"Saved {outroot}.reverse.png (empty - no non-zero entries)")
            else:
                ys_r, xs_r = np.nonzero(nz_mask)
                vals_r = arr[ys_r, xs_r].astype(np.float64)
                if args.logscale:
                    vals_r = np.log1p(vals_r)
                if args.vmax is None:
                    try:
                        auto_vmax = float(np.percentile(vals_r, 99))
                        if auto_vmax <= 0:
                            auto_vmax = float(vals_r.max())
                        if auto_vmax <= 0:
                            auto_vmax = 1.0
                    except Exception:
                        auto_vmax = float(vals_r.max()) if vals_r.size else 1.0
                    vmax_use = auto_vmax
                else:
                    vmax_use = args.vmax
                cmap_obj = plt.get_cmap(args.cmap)
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
                if args.size_by_value:
                    vmax_val = vals_r.max() if vals_r.size else 1.0
                    rel = (vals_r / vmax_val)
                    sizes = (np.sqrt(rel) * args.marker_size * args.size_scale) + (args.marker_size * 0.1)
                else:
                    sizes = np.full_like(vals_r, fill_value=args.marker_size, dtype=float)
                fig, ax = plt.subplots(figsize=(min(max(4, ncols/200*10),20), min(max(4, nrows/200*10),20)))
                sc = ax.scatter(xs_r, ys_r, c=vals_r, cmap=cmap_with_transparency, s=sizes, vmax=vmax_use,
                                alpha=args.alpha, marker='s', edgecolors='none')
                cbar = fig.colorbar(sc, ax=ax)
                cbar.set_label('log1p(counts)' if args.logscale else 'counts')
                ax.set_title(f"Reverse matches (k={args.k}) - {root_file_name_sans_dir(args.reference)}")
                def choose_ticks(n, max_ticks):
                    if n <= max_ticks:
                        return np.arange(n)
                    locs = np.linspace(0, n - 1, num=max_ticks, dtype=int)
                    return np.unique(locs)
                xticks = choose_ticks(ncols, 10)
                yticks = choose_ticks(nrows, 10)
                if x_coords is not None:
                    x_coords_arr = np.asarray(x_coords)
                    if x_coords_arr.shape[0] == ncols:
                        xtick_labels = [str(int(x_coords_arr[x])) for x in xticks]
                        ax.set_xlabel('reference base (k-mer start, 1-based)')
                    else:
                        xtick_labels = [str(int(x) + 1) for x in xticks]
                        ax.set_xlabel('reference k-mer start (1-based)')
                else:
                    xtick_labels = [str(int(x) + 1) for x in xticks]
                    ax.set_xlabel('reference k-mer start (1-based)')
                ytick_labels = [str(int(y) + 1) for y in yticks]
                ax.set_xticks(xticks)
                ax.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
                ax.set_yticks(yticks)
                ax.set_yticklabels(ytick_labels, fontsize=8)
                ax.set_ylabel(f"comp (bin size={bin_size_comp}) k-mer start (1-based)")
                plt.tight_layout()
                plt.savefig(f"{outroot}.reverse.png", dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"Saved {outroot}.reverse.png")

    print("Done.")

if __name__ == "__main__":
    main()
