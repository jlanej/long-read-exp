#!/usr/bin/env python3
"""
Position-agnostic, aligner-free deletion discovery via conserved k-mer anchors.

- Input: multi-fasta (all sequences ~30 kb, same orientation)
- Steps:
  1) Index all k-mers (k=31 by default) for each sequence; record first position(s).
  2) Select "anchors": k-mers present in >= anchor_presence fraction of samples AND
     single-copy in >= single_copy_fraction of those samples.
  3) Order anchors by their median position across samples that contain them.
  4) For adjacent ordered anchors (and also a stride>1 option), compute inter-anchor distance per sample.
  5) Score each anchor pair for bimodality: check variance & a simple dip-like metric;
     call a "short" cluster vs "long" cluster via 1D k-means (k=2).
  6) Pick the best-discriminating pair; call per-sample deletion if assigned to "short" cluster.
  7) Plot distance heatmap (samples x anchor-pairs) with enhanced visualizations.

Outputs:
  - anchors.csv
  - pair_calls.csv
  - sample_calls.csv
  - distance_heatmap.png (basic plot)
  - enhanced_heatmap.png (multi-panel diagnostic plot with deletion highlighting)
  - summary_plots.png (statistical overview plots)

Enhanced plotting features:
  - Multi-panel heatmap showing raw distances, z-scores, deletion calls, and distributions
  - Color-coded deletion events (red) vs normal samples (green)
  - Statistical overlays showing cluster centers and separation metrics
  - Summary plots with pair quality metrics and deletion call statistics
  - Backward compatible with original basic plotting mode
"""

import argparse, sys, math
from collections import defaultdict, Counter
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.vq import kmeans2
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import zscore

def rc(seq):
    tbl = str.maketrans("ACGTNacgtn","TGCANtgcan")
    return seq.translate(tbl)[::-1]

def all_kmer_positions(seq, k):
    """Return dict kmer -> first position (and flag if multi-copy)."""
    pos = {}
    counts = Counter()
    L = len(seq)
    for i in range(L - k + 1):
        kmer = seq[i:i+k]
        counts[kmer] += 1
        if kmer not in pos:
            pos[kmer] = i
    multi = {k for k,c in counts.items() if c>1}
    return pos, multi

def load_sequences(fasta):
    seqs = []
    for rec in SeqIO.parse(fasta, "fasta"):
        s = str(rec.seq).upper()
        seqs.append((rec.id, s))
    return seqs

def choose_anchors(seqs, k=31, presence=0.90, single_copy_fraction=0.95, max_anchors=3000):
    """
    Find kmers that are present in >=presence fraction of sequences AND
    single-copy in >=single_copy_fraction of sequences where present.
    """
    n = len(seqs)
    kmer_presence = Counter()
    kmer_singlecopy_ok = Counter()
    kmer_firstpos = defaultdict(list)

    # pass 1: per-seq index
    perseq_pos = {}
    perseq_multi = {}
    for sid, s in seqs:
        pos, multi = all_kmer_positions(s, k)
        perseq_pos[sid] = pos
        perseq_multi[sid] = multi
        for kmer in pos.keys():
            kmer_presence[kmer] += 1
            if kmer not in multi:
                kmer_singlecopy_ok[kmer] += 1

    # filter anchors
    min_present = math.ceil(presence * n)
    anchors = []
    for kmer, cnt in kmer_presence.items():
        if cnt >= min_present:
            sc_ok = kmer_singlecopy_ok[kmer]
            frac_sc = sc_ok / cnt
            if frac_sc >= single_copy_fraction:
                anchors.append(kmer)

    # cap and return structures
    if len(anchors) > max_anchors:
        # keep the most prevalent anchors (ties arbitrary), favoring those with highest single-copy fraction
        anchors = sorted(anchors, key=lambda km: (kmer_presence[km], kmer_singlecopy_ok[km]/kmer_presence[km]), reverse=True)[:max_anchors]

    # build position table for anchors
    rows = []
    for km in anchors:
        rows.append({
            "kmer": km,
            "present_in": kmer_presence[km],
            "single_copy_ok": kmer_singlecopy_ok[km],
            "present_frac": kmer_presence[km]/n,
            "single_copy_frac": kmer_singlecopy_ok[km]/kmer_presence[km]
        })
    
    if len(rows) == 0:
        print(f"No anchors found with presence >= {presence} and single_copy_frac >= {single_copy_fraction}", file=sys.stderr)
        return pd.DataFrame(columns=["kmer","present_in","single_copy_ok","present_frac","single_copy_frac"]), perseq_pos, perseq_multi
    
    anchors_df = pd.DataFrame(rows).sort_values(["present_in","single_copy_frac"], ascending=[False,False])
    return anchors_df, perseq_pos, perseq_multi

def order_anchors_by_median_pos(anchors_df, seqs, perseq_pos):
    # compute median position per anchor across sequences where present
    med_pos = []
    for km in anchors_df["kmer"]:
        ps = []
        for sid,_ in seqs:
            p = perseq_pos[sid].get(km, None)
            if p is not None:
                ps.append(p)
        if len(ps) == 0:
            med = np.nan
        else:
            med = float(np.median(ps))
        med_pos.append(med)
    anchors_df = anchors_df.copy()
    anchors_df["median_pos"] = med_pos
    anchors_df = anchors_df.dropna().sort_values("median_pos").reset_index(drop=True)
    return anchors_df

def compute_distances(seqs, anchors_df, perseq_pos, pair_stride=1):
    """
    For ordered anchors, compute distances between anchor i and i+stride for each sequence.
    Returns a long-format dataframe: sequence_id, pair_index, anchor_i, anchor_j, distance (or NaN if missing).
    """
    ordered = anchors_df["kmer"].tolist()
    pairs = []
    for i in range(len(ordered) - pair_stride):
        a = ordered[i]
        b = ordered[i + pair_stride]
        pairs.append((i, a, b))
    rows = []
    for sid,_ in seqs:
        pmap = perseq_pos[sid]
        for idx, a, b in pairs:
            pa = pmap.get(a, None)
            pb = pmap.get(b, None)
            if pa is None or pb is None:
                dist = np.nan
            else:
                dist = pb - pa
            rows.append((sid, idx, a, b, dist))
    df = pd.DataFrame(rows, columns=["id","pair_idx","anchor_a","anchor_b","distance"])
    return df

def score_bimodality_and_call(df_pair):
    """
    Given distances for one anchor pair across samples, check for bimodality and split by k-means (k=2).
    Return dict with stats; add per-sample label ('short'/'long'/NaN).
    """
    x = df_pair["distance"].dropna().values.astype(float)
    if len(x) < 10 or np.nanstd(x) < 50:  # too few or too tight to be interesting
        return None, None
    # k-means k=2 on 1D
    cents, labels = kmeans2(x.reshape(-1,1), 2, minit='points')
    # label centers as short/long
    short_center = np.min(cents)
    long_center  = np.max(cents)
    # soft rule: require good separation
    sep = long_center - short_center
    if sep < 500:  # require at least ~500 bp separation to consider (tune as needed)
        return None, None
    # assign labels back
    arr_labels = np.array(["short" if abs(v-short_center) < abs(v-long_center) else "long" for v in x])
    # basic quality metrics
    frac_short = (arr_labels=="short").mean()
    frac_long = 1 - frac_short
    stats = {
        "short_center": float(short_center),
        "long_center": float(long_center),
        "separation": float(sep),
        "frac_short": float(frac_short),
        "frac_long": float(frac_long),
        "n": int(len(x)),
        "std": float(np.std(x))
    }
    # build per-sample mapping
    mapping = dict(zip(df_pair["id"].dropna(), arr_labels))
    return stats, mapping

def create_enhanced_heatmap(M, row_ids_sorted, col_labels, calls_df, pair_df, top_pairs, out_prefix):
    """
    Create enhanced heatmaps with better visualization and deletion highlighting.
    """
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("viridis")
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Enhanced K-mer Anchor Distance Analysis', fontsize=16, fontweight='bold')
    
    # 1. Enhanced distance heatmap with better color scheme
    ax1 = axes[0, 0]
    
    # Create mask for NaN values
    mask = np.isnan(M)
    
    # Use seaborn heatmap for better visualization
    sns.heatmap(M, 
                mask=mask,
                cmap='viridis',
                cbar_kws={'label': 'Inter-anchor distance (bp)'},
                xticklabels=[f"Pair {i}" for i in col_labels],
                yticklabels=False,  # Too many samples to show labels
                ax=ax1,
                robust=True)
    
    ax1.set_title('Distance Heatmap (Raw Values)', fontweight='bold')
    ax1.set_xlabel('Anchor Pairs')
    ax1.set_ylabel('Samples (sorted by deletion call)')
    
    # 2. Z-score normalized heatmap to highlight patterns
    ax2 = axes[0, 1]
    
    # Calculate z-scores column-wise (per anchor pair)
    M_zscore = np.full_like(M, np.nan)
    for j in range(M.shape[1]):
        col_data = M[:, j]
        valid_data = col_data[~np.isnan(col_data)]
        if len(valid_data) > 1:
            z_scores = zscore(valid_data)
            M_zscore[~np.isnan(col_data), j] = z_scores
    
    mask_z = np.isnan(M_zscore)
    sns.heatmap(M_zscore,
                mask=mask_z,
                cmap='RdBu_r',
                center=0,
                cbar_kws={'label': 'Z-score'},
                xticklabels=[f"Pair {i}" for i in col_labels],
                yticklabels=False,
                ax=ax2,
                vmin=-3, vmax=3)
    
    ax2.set_title('Z-score Normalized Distances', fontweight='bold')
    ax2.set_xlabel('Anchor Pairs')
    ax2.set_ylabel('Samples')
    
    # 3. Binary deletion calls heatmap
    ax3 = axes[1, 0]
    
    # Create binary matrix based on calls
    call_order = [row_ids_sorted[i] for i in range(len(row_ids_sorted))]
    binary_matrix = np.full((len(call_order), len(top_pairs)), np.nan)
    
    for i, sample_id in enumerate(call_order):
        sample_call = calls_df[calls_df['id'] == sample_id]['call'].values
        if len(sample_call) > 0:
            call = sample_call[0]
            for j in range(len(top_pairs)):
                if call == 'deletion':
                    binary_matrix[i, j] = 1
                elif call == 'no_deletion':
                    binary_matrix[i, j] = 0
                else:  # unknown
                    binary_matrix[i, j] = 0.5
    
    # Custom colormap for binary calls
    colors = ['#2E8B57', '#FFE4B5', '#DC143C']  # no_deletion, unknown, deletion
    custom_cmap = sns.blend_palette(colors, n_colors=100, as_cmap=True)
    
    mask_binary = np.isnan(binary_matrix)
    sns.heatmap(binary_matrix,
                mask=mask_binary,
                cmap=custom_cmap,
                cbar_kws={'label': 'Deletion Call', 
                         'ticks': [0, 0.5, 1],
                         'format': plt.FuncFormatter(lambda x, p: {0: 'No Deletion', 0.5: 'Unknown', 1: 'Deletion'}[x])},
                xticklabels=[f"Pair {i}" for i in col_labels],
                yticklabels=False,
                ax=ax3,
                vmin=0, vmax=1)
    
    ax3.set_title('Deletion Calls by Sample', fontweight='bold')
    ax3.set_xlabel('Anchor Pairs')
    ax3.set_ylabel('Samples')
    
    # 4. Distance distribution plot highlighting deletions
    ax4 = axes[1, 1]
    
    # Get the best discriminating pair (first one)
    if len(top_pairs) > 0:
        best_pair_idx = top_pairs[0]
        best_pair_distances = M[:, 0]  # First column corresponds to best pair
        
        # Separate deletion vs no-deletion samples
        deletion_distances = []
        no_deletion_distances = []
        
        for i, sample_id in enumerate(call_order):
            sample_call = calls_df[calls_df['id'] == sample_id]['call'].values
            if len(sample_call) > 0 and not np.isnan(best_pair_distances[i]):
                call = sample_call[0]
                if call == 'deletion':
                    deletion_distances.append(best_pair_distances[i])
                elif call == 'no_deletion':
                    no_deletion_distances.append(best_pair_distances[i])
        
        # Create histogram
        bins = np.linspace(np.nanmin(best_pair_distances), np.nanmax(best_pair_distances), 30)
        
        ax4.hist(no_deletion_distances, bins=bins, alpha=0.7, label='No Deletion', 
                color='#2E8B57', density=True)
        ax4.hist(deletion_distances, bins=bins, alpha=0.7, label='Deletion', 
                color='#DC143C', density=True)
        
        ax4.set_xlabel('Distance (bp)')
        ax4.set_ylabel('Density')
        ax4.set_title(f'Distance Distribution (Best Pair {best_pair_idx})', fontweight='bold')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        # Add vertical lines for cluster centers if available
        best_stats = pair_df[pair_df['pair_idx'] == best_pair_idx]
        if len(best_stats) > 0:
            short_center = best_stats['short_center'].iloc[0]
            long_center = best_stats['long_center'].iloc[0]
            ax4.axvline(short_center, color='red', linestyle='--', alpha=0.8, 
                       label=f'Short Center ({short_center:.0f} bp)')
            ax4.axvline(long_center, color='green', linestyle='--', alpha=0.8,
                       label=f'Long Center ({long_center:.0f} bp)')
            ax4.legend()
    
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_enhanced_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Wrote {out_prefix}_enhanced_heatmap.png", file=sys.stderr)

def create_summary_plots(pair_df, calls_df, out_prefix):
    """
    Create summary plots showing overall statistics and patterns.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('K-mer Anchor Analysis Summary', fontsize=16, fontweight='bold')
    
    # 1. Pair separation vs standard deviation
    ax1 = axes[0, 0]
    scatter = ax1.scatter(pair_df['separation'], pair_df['std'], 
                         c=pair_df['frac_short'], cmap='RdYlBu_r',
                         alpha=0.7, s=60)
    ax1.set_xlabel('Separation (bp)')
    ax1.set_ylabel('Standard Deviation')
    ax1.set_title('Pair Quality Metrics')
    ax1.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax1, label='Fraction with Deletions')
    
    # 2. Distribution of separation distances
    ax2 = axes[0, 1]
    ax2.hist(pair_df['separation'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax2.set_xlabel('Separation Distance (bp)')
    ax2.set_ylabel('Number of Pairs')
    ax2.set_title('Distribution of Pair Separations')
    ax2.grid(True, alpha=0.3)
    
    # 3. Deletion call summary
    ax3 = axes[1, 0]
    call_counts = calls_df['call'].value_counts()
    colors = {'deletion': '#DC143C', 'no_deletion': '#2E8B57', 'unknown': '#FFE4B5'}
    pie_colors = [colors.get(call, 'gray') for call in call_counts.index]
    
    ax3.pie(call_counts.values, labels=call_counts.index, autopct='%1.1f%%',
           colors=pie_colors, startangle=90)
    ax3.set_title('Deletion Call Summary')
    
    # 4. Top pairs ranking
    ax4 = axes[1, 1]
    top_10 = pair_df.head(10).copy()
    top_10['rank'] = range(1, len(top_10) + 1)
    
    bars = ax4.bar(top_10['rank'], top_10['separation'], 
                   color=plt.cm.viridis(top_10['frac_short']))
    ax4.set_xlabel('Pair Rank')
    ax4.set_ylabel('Separation (bp)')
    ax4.set_title('Top 10 Discriminative Pairs')
    ax4.set_xticks(top_10['rank'])
    ax4.grid(True, alpha=0.3)
    
    # Add text annotations for separation values
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_summary_plots.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Wrote {out_prefix}_summary_plots.png", file=sys.stderr)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seqs", required=True, help="multi-fasta of ~250 sequences")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--presence", type=float, default=0.90, help="anchor must be present in >= this fraction of samples")
    ap.add_argument("--single-copy-frac", type=float, default=0.95, help="anchor must be single-copy in >= this fraction of present samples")
    ap.add_argument("--max-anchors", type=int, default=3000)
    ap.add_argument("--pair-stride", type=int, default=1, help="compute distances between anchor i and i+stride")
    ap.add_argument("--heatmap-top", type=int, default=50, help="plot top-N discriminative pairs")
    ap.add_argument("--out-prefix", default="pos_agnostic")
    ap.add_argument("--enhanced-plots", action="store_true", help="create enhanced diagnostic plots with better visualizations")
    ap.add_argument("--basic-plots-only", action="store_true", help="create only the basic distance heatmap (legacy mode)")
    args = ap.parse_args()

    seqs = load_sequences(args.seqs)
    print(f"Loaded {len(seqs)} sequences.", file=sys.stderr)

    anchors_df, perseq_pos, perseq_multi = choose_anchors(
        seqs, k=args.k, presence=args.presence, single_copy_fraction=args.single_copy_frac, max_anchors=args.max_anchors
    )
    anchors_df.to_csv(args.out_prefix + "_anchors.csv", index=False)
    print(f"Anchors selected: {len(anchors_df)} (wrote anchors.csv)", file=sys.stderr)

    ordered_df = order_anchors_by_median_pos(anchors_df, seqs, perseq_pos)
    print(f"Ordered anchors retained: {len(ordered_df)}", file=sys.stderr)

    dist_df = compute_distances(seqs, ordered_df, perseq_pos, pair_stride=args.pair_stride)

    # score each pair
    pair_stats = []
    per_pair_labels = {}
    for idx, sub in dist_df.groupby("pair_idx"):
        stats, labels = score_bimodality_and_call(sub)
        if stats is None:
            continue
        a = sub["anchor_a"].iloc[0]; b = sub["anchor_b"].iloc[0]
        stats.update({"pair_idx": int(idx), "anchor_a": a, "anchor_b": b})
        pair_stats.append(stats)
        per_pair_labels[idx] = labels

    if not pair_stats:
        print("No strongly bimodal anchor pairs found. Consider lowering k, relaxing presence, or increasing max_anchors.", file=sys.stderr)
        pair_df = pd.DataFrame(columns=["pair_idx","anchor_a","anchor_b","short_center","long_center","separation","frac_short","frac_long","n","std"])
        pair_df.to_csv(args.out_prefix + "_pair_calls.csv", index=False)
        sys.exit(0)

    pair_df = pd.DataFrame(pair_stats).sort_values(["separation","std","n"], ascending=[False,False,False])
    pair_df.to_csv(args.out_prefix + "_pair_calls.csv", index=False)
    print(f"Scored {len(pair_df)} bimodal pairs (wrote pair_calls.csv)", file=sys.stderr)

    # choose best pair for per-sample final call (largest separation)
    best = pair_df.iloc[0]
    best_idx = int(best["pair_idx"])
    best_labels = per_pair_labels[best_idx]
    # assemble per-sample calls
    calls = []
    for sid,_ in seqs:
        lab = best_labels.get(sid, None)
        calls.append({"id": sid, "call": ("deletion" if lab=="short" else ("no_deletion" if lab=="long" else "unknown")),
                      "evidence_pair_idx": best_idx})
    calls_df = pd.DataFrame(calls)
    calls_df.to_csv(args.out_prefix + "_sample_calls.csv", index=False)
    print(f"Wrote per-sample calls to sample_calls.csv using pair_idx={best_idx}", file=sys.stderr)

    # Heatmap: pick top-N pairs and plot distances
    top_pairs = pair_df.head(args.heatmap_top)["pair_idx"].tolist()
    mat_rows = []
    row_ids = [sid for sid,_ in seqs]
    col_labels = []
    for pidx in top_pairs:
        sub = dist_df[dist_df["pair_idx"]==pidx].set_index("id")
        col = []
        for sid in row_ids:
            v = sub.loc[sid]["distance"] if sid in sub.index else np.nan
            col.append(v)
        mat_rows.append(col)
        a = dist_df[dist_df["pair_idx"]==pidx]["anchor_a"].iloc[0]
        b = dist_df[dist_df["pair_idx"]==pidx]["anchor_b"].iloc[0]
        col_labels.append(f"{pidx}")

    M = np.array(mat_rows).T  # samples x pairs
    # simple row order: group by best call
    order = np.argsort([{"deletion":0,"unknown":1,"no_deletion":2}[x] for x in calls_df["call"]])
    M = M[order,:]
    row_ids_sorted = [row_ids[i] for i in order]

    # Create plots based on command line options
    if not args.basic_plots_only:
        # Create enhanced diagnostic plots by default
        if len(pair_stats) > 0:
            create_enhanced_heatmap(M, row_ids_sorted, col_labels, calls_df, pair_df, top_pairs, args.out_prefix)
            create_summary_plots(pair_df, calls_df, args.out_prefix)
        else:
            print("No bimodal pairs found - creating basic plots only", file=sys.stderr)
        
        if args.enhanced_plots:
            print("Enhanced plots already created as default", file=sys.stderr)

    # Always create the basic plot for backward compatibility, unless enhanced-only mode
    if args.basic_plots_only or not args.enhanced_plots or len(pair_stats) == 0:
        plt.figure(figsize=(max(8, len(top_pairs)*0.25), max(6, len(row_ids_sorted)*0.02)))
        plt.imshow(M, aspect='auto', interpolation='nearest')
        plt.colorbar(label="Inter-anchor distance (bp)")
        plt.xlabel("Top bimodal anchor pairs (index)")
        plt.ylabel("Samples (sorted by call)")
        plt.title("Aligner-free, position-agnostic deletion signal via k-mer anchors")
        plt.tight_layout()
        plt.savefig(args.out_prefix + "_distance_heatmap.png", dpi=200)
        print("Wrote distance_heatmap.png", file=sys.stderr)

if __name__ == "__main__":
    main()
