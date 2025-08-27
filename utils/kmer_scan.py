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
import csv


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
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()

    for ax, pidx in zip(axes, top_pairs):
        sub = dist_df[dist_df["pair_idx"] == pidx].dropna(subset=["distance"])
        sub = sub.merge(calls_df[["id", "call"]], on="id", how="left")

        # Plot histogram of all distances
        ax.hist(sub["distance"], bins=50, alpha=0.6, color="gray", label="all samples")

        # Overlay deletion distances
        del_dists = sub[sub["call"] == "deletion"]["distance"].values
        if len(del_dists) > 0:
            for d in del_dists:
                ax.axvline(d, color="red", linestyle="--", linewidth=1.2, alpha=0.8)

        # Use anchor names from dist_df
        a = dist_df.loc[dist_df["pair_idx"] == pidx, "anchor_a"].iloc[0]
        b = dist_df.loc[dist_df["pair_idx"] == pidx, "anchor_b"].iloc[0]

        ax.set_title(f"Pair {pidx}: {a[:6]}… → {b[:6]}…")
        ax.set_xlabel("Distance (bp)")
        ax.set_ylabel("Count")

    # Hide unused axes
    for ax in axes[len(top_pairs):]:
        ax.axis("off")

    fig.suptitle(f"Distance distributions for top {len(top_pairs)} bimodal anchor pairs", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(args.out_prefix + "_distance_histograms.png", dpi=200)
    print("Wrote distance_histograms.png", file=sys.stderr)

    # -----------------------
    # Scatter + Heatmap of anchor-pair distance variability (no npz saving)
    # -----------------------

    # Build mapping of anchor -> median_pos (from ordered_df)
    anchor_to_pos = dict(zip(ordered_df["kmer"], ordered_df["median_pos"]))

    # 1) Scatterplot: one point per pair_idx present in dist_df
    scatter_rows = []
    for idx, sub in dist_df.groupby("pair_idx"):
        dists = sub["distance"].dropna().values.astype(float)
        if len(dists) < 2:
            continue
        stdv = float(np.std(dists))
        a = sub["anchor_a"].iloc[0]
        b = sub["anchor_b"].iloc[0]
        if a not in anchor_to_pos or b not in anchor_to_pos:
            continue
        pos_a = anchor_to_pos[a]
        pos_b = anchor_to_pos[b]
        scatter_rows.append((pos_a, pos_b, stdv, idx, a, b))

    if len(scatter_rows) == 0:
        print("No anchor pairs with enough data to plot scatter/heatmap.", file=sys.stderr)
    else:
        X = np.array([r[0] for r in scatter_rows])
        Y = np.array([r[1] for r in scatter_rows])
        C = np.array([r[2] for r in scatter_rows])

        plt.figure(figsize=(8,8))
        sc = plt.scatter(X, Y, c=C, cmap="viridis", s=45, edgecolors="k", linewidths=0.2)
        cb = plt.colorbar(sc)
        cb.set_label("Std. dev. of inter-anchor distance (bp)")
        plt.xlabel("Anchor A median position (bp)")
        plt.ylabel("Anchor B median position (bp)")
        plt.title("Anchor-pair variability (std of distances)")
        plt.grid(alpha=0.2)
        plt.tight_layout()
        out_scatter = args.out_prefix + "_anchorpair_variability_scatter.png"
        plt.savefig(out_scatter, dpi=300)
        plt.close()
        print(f"Wrote {out_scatter}", file=sys.stderr)

        # -----------------------
        # Junction-finding & breakpoint plotting for deletion-carrying samples
        # Place after pair_df, per_pair_labels, seqs, perseq_pos, ordered_df exist
        # -----------------------

        def find_junction_in_sample(sample_seq, pa, pb, anchor_k, k_j=21, max_left=200, max_right=200):
            """
            Attempt to find an exact junction (left_suffix + right_prefix) in sample_seq.
            - pa = position of anchor A (start of anchor kmer) in sample_seq
            - pb = position of anchor B (start of anchor kmer) in sample_seq
            - anchor_k = anchor k-mer length (args.k)
            - k_j = junction k-mer length to try (<= anchor_k recommended)
            Returns:
              (breakpoint_pos, match_pos, junction_kmer, left_len, right_len)
              where breakpoint_pos is the index in sample_seq where left/right join occurs (position of first base of right part)
              or None if no exact junction found.
            """
            s = sample_seq
            L = len(s)

            # Build left context ending at pa + anchor_k (end of anchor A)
            left_end = min(L, pa + anchor_k)
            left_start = max(0, left_end - max_left - k_j)  # keep enough left context
            left_ctx = s[left_start:left_end]

            # Build right context starting at pb (start of anchor B)
            right_start = max(0, pb)
            right_end = min(L, pb + max_right + k_j)
            right_ctx = s[right_start:right_end]

            # boundary for search: whole sample sequence is fine (we want absolute match)
            # Try splits: left_take from 1..k_j (left tail length)
            for left_take in range(k_j, 0, -1):  # prefer longer left tail first
                right_take = k_j - left_take
                if right_take < 0:
                    continue
                left_part = left_ctx[-left_take:] if left_take <= len(left_ctx) else None
                right_part = right_ctx[:right_take] if right_take == 0 or right_take <= len(right_ctx) else None
                if left_part is None or right_part is None:
                    continue
                junction = left_part + right_part
                # search for junction in sample sequence
                pos = s.find(junction)
                if pos != -1:
                    # breakpoint is pos + len(left_part)  (first base of the right part)
                    bp = pos + len(left_part)
                    return int(bp), int(pos), junction, int(left_take), int(right_take)
            return None, None, None, None, None

        # helper: map seq id -> seq string
        seq_dict = {sid: s for sid, s in seqs}
        anchor_to_medpos = dict(zip(ordered_df["kmer"], ordered_df["median_pos"]))

        TOP_K = min(20, len(pair_df))  # adjust: how many top pairs to inspect
        out_dir_prefix = args.out_prefix

        # ensure output dir exists
        import os
        out_dir = out_dir_prefix + "_breakpoints"
        os.makedirs(out_dir, exist_ok=True)

        # iterate top pairs
        for _, prow in pair_df.head(TOP_K).iterrows():
            pidx = int(prow["pair_idx"])
            a_kmer = prow["anchor_a"]
            b_kmer = prow["anchor_b"]
            a_med = anchor_to_medpos.get(a_kmer, np.nan)
            b_med = anchor_to_medpos.get(b_kmer, np.nan)

            # per-sample labels for this pair
            labels = per_pair_labels.get(pidx, {})
            # collect results
            rows = []
            for sid, s in seqs:
                lab = labels.get(sid, None)
                if lab != "short":  # only analyze deletion-carrying samples (short)
                    continue
                pmap = perseq_pos.get(sid, {})
                pa = pmap.get(a_kmer, None)
                pb = pmap.get(b_kmer, None)
                if pa is None or pb is None:
                    # anchor missing — we cannot localize; skip or note as missing
                    rows.append({
                        "sample": sid, "status": "missing_anchor", "pa": pa, "pb": pb,
                        "breakpoint_sample_bp": np.nan, "match_pos": np.nan, "junction": "", "left_take": np.nan,
                        "right_take": np.nan
                    })
                    continue

                # try to find junction exact sequence using k_j (try multiple lengths)
                found = False
                for k_j in (min(args.k, 31), 25, 21, 17):  # try decreasing k_j to be permissive
                    bp_pos, match_pos, junction, left_take, right_take = find_junction_in_sample(
                        seq_dict[sid], pa, pb, anchor_k=args.k, k_j=k_j, max_left=400, max_right=400
                    )
                    if bp_pos is not None:
                        rows.append({
                            "sample": sid,
                            "status": "junction_found",
                            "pa": pa, "pb": pb,
                            "breakpoint_sample_bp": bp_pos,
                            "match_pos": match_pos,
                            "junction": junction,
                            "k_j": k_j,
                            "left_take": left_take,
                            "right_take": right_take
                        })
                        found = True
                        break
                if not found:
                    # fallback: report ambiguous interval between pa+anchor_k and pb-1 in sample coordinates
                    left_bound = pa + args.k
                    right_bound = pb - 1
                    if right_bound < left_bound:
                        # if pb <= pa+anchor_k, set narrow ambiguous
                        left_bound = max(0, pa + args.k - 5)
                        right_bound = min(len(seq_dict[sid]) - 1, pb + 5)
                    rows.append({
                        "sample": sid,
                        "status": "no_junction_found",
                        "pa": pa, "pb": pb,
                        "breakpoint_sample_bp": np.nan,
                        "ambig_start": left_bound,
                        "ambig_end": right_bound,
                        "junction": "",
                        "k_j": np.nan
                    })

            # write CSV summarizing per-sample breakpoint findings for this pair
            csv_path = os.path.join(out_dir, f"pair_{pidx}_breakpoints.csv")
            # gather all fieldnames dynamically
            all_keys = set()
            for r in rows:
                all_keys.update(r.keys())
            fieldnames = sorted(list(all_keys))
            with open(csv_path, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)
            print("Wrote", csv_path, file=sys.stderr)

            # build plots for this pair
            # 1) histogram of found breakpoints (sample bp)
            found_bps = [r["breakpoint_sample_bp"] for r in rows if r.get("status") == "junction_found"]
            ambigs = [(r.get("ambig_start"), r.get("ambig_end")) for r in rows if
                      r.get("status") == "no_junction_found"]
            missing = [r["sample"] for r in rows if r.get("status") == "missing_anchor"]

            plt.figure(figsize=(8, 4))
            if len(found_bps) > 0:
                plt.hist(found_bps, bins=30, alpha=0.7, color='C2', label='exact junction bp')
            if len(ambigs) > 0:
                # draw ambiguous intervals as horizontal tick bars
                for (s, e) in ambigs:
                    if np.isnan(s) or np.isnan(e):
                        continue
                    plt.axvspan(s, e, color='orange', alpha=0.15)
            # overlay pa+args.k and pb for one representative sample (use median of pa/pb across deletion samples)
            pa_vals = [r["pa"] for r in rows if r.get("pa") is not None]
            pb_vals = [r["pb"] for r in rows if r.get("pb") is not None]
            if len(pa_vals) > 0:
                plt.axvline(np.median(pa_vals) + args.k, color='k', linestyle='--', label='median pa+anchor_k')
            if len(pb_vals) > 0:
                plt.axvline(np.median(pb_vals), color='gray', linestyle=':', label='median pb')
            plt.xlabel("Breakpoint coordinate (bp) in sample sequence")
            plt.ylabel("Count / coverage")
            plt.title(
                f"Pair {pidx} junction positions across deletion samples\nanchors medpos [{a_med:.0f} -> {b_med:.0f}]")
            plt.legend()
            out_hist = os.path.join(out_dir, f"pair_{pidx}_breakpoint_hist.png")
            plt.tight_layout()
            plt.savefig(out_hist, dpi=200)
            plt.close()
            print("Wrote", out_hist, file=sys.stderr)

            # 2) per-sample scatter (y = sample index; x = breakpoint bp or ambiguous interval midpoint)
            sample_positions = []
            sample_labels = []
            for r in rows:
                sid = r["sample"]
                if r["status"] == "junction_found":
                    sample_positions.append(r["breakpoint_sample_bp"])
                elif r["status"] == "no_junction_found":
                    s, e = r["ambig_start"], r["ambig_end"]
                    sample_positions.append((s + e) / 2.0 if (not np.isnan(s) and not np.isnan(e)) else np.nan)
                else:
                    sample_positions.append(np.nan)
                sample_labels.append(sid)
            y = np.arange(len(sample_labels))
            plt.figure(figsize=(8, max(4, len(sample_labels) * 0.06)))
            plt.scatter(sample_positions, y, c='C1', s=16)
            # annotate missing anchors
            for i, r in enumerate(rows):
                if r["status"] == "missing_anchor":
                    plt.scatter(-10, i, marker='x', color='red')  # off-left marker for missing
            plt.yticks(y, sample_labels, fontsize=6)
            plt.xlabel("Breakpoint coordinate (bp) in sample")
            plt.ylabel("Sample")
            plt.title(f"Per-sample breakpoint positions for pair {pidx}")
            plt.tight_layout()
            out_scatter = os.path.join(out_dir, f"pair_{pidx}_breakpoint_scatter.png")
            plt.savefig(out_scatter, dpi=200)
            plt.close()
            print("Wrote", out_scatter, file=sys.stderr)

        # done for top pairs
        print("Junction search completed. CSVs and PNGs are in", out_dir, file=sys.stderr)


if __name__ == "__main__":
    main()
