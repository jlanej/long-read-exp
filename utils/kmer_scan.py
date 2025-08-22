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
  7) Plot distance heatmap (samples x anchor-pairs).

Outputs:
  - anchors.csv
  - pair_calls.csv
  - sample_calls.csv
  - distance_heatmap.png
"""

import argparse, sys, math
from collections import defaultdict, Counter
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2

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

    plt.figure(figsize=(max(8, len(top_pairs)*0.25), max(6, len(row_ids_sorted)*0.02)))
    plt.imshow(M, aspect='auto', interpolation='nearest')
    plt.colorbar(label="Inter-anchor distance (bp)")
    plt.xlabel("Top bimodal anchor pairs (index)")
    plt.ylabel("Samples (sorted by call)")
    plt.title("Aligner-free, position-agnostic deletion signal via k-mer anchors")
    plt.tight_layout()
    plt.savefig(args.out_prefix + "_distance_heatmap.png", dpi=200)
    print("Wrote distance_heatmap.png", file=sys.stderr)

    # Histograms for the top-N bimodal pairs
    top_pairs = pair_df.head(args.heatmap_top)["pair_idx"].tolist()

    ncols = 3
    nrows = int(math.ceil(len(top_pairs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes = axes.flatten()

    for ax, pidx in zip(axes, top_pairs):
        sub = dist_df[dist_df["pair_idx"] == pidx].dropna(subset=["distance"])
        sub = sub.merge(calls_df[["id", "call"]], on="id", how="left")

        # Plot histogram of all distances
        ax.hist(sub["distance"], bins=50, alpha=0.6, color="gray", label="all samples")

        # Overlay deletion distances as vertical red lines
        del_dists = sub[sub["call"]=="deletion"]["distance"].values
        for d in del_dists:
            ax.axvline(d, color="red", linestyle="--", linewidth=1.2, alpha=0.8)

        a = sub["anchor_a"].iloc[0]
        b = sub["anchor_b"].iloc[0]
        ax.set_title(f"Pair {pidx}\n{a[:6]}… → {b[:6]}…")
        ax.set_xlabel("Distance (bp)")
        ax.set_ylabel("Count")

    # Remove unused axes
    for ax in axes[len(top_pairs):]:
        ax.axis("off")

    fig.suptitle(f"Distance distributions for top {len(top_pairs)} bimodal anchor pairs", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(args.out_prefix + "_distance_histograms.png", dpi=200)
    print("Wrote distance_histograms.png", file=sys.stderr)

if __name__ == "__main__":
    main()
