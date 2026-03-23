#!/usr/bin/env python3
"""
Merges per-sample diversity outputs and computes alpha and beta diversity
across all samples simultaneously.

Optional --w_seg accepts a comma-separated list of segment tokens, each
of which may contain pipe-separated aliases that are treated as a single
group. Example:
    --w_seg "AWF,CL,LSW|LSF"
produces Whittaker beta estimates for:
    - all samples (global)
    - samples whose sample_id contains "AWF"
    - samples whose sample_id contains "CL"
    - samples whose sample_id contains "LSW" OR "LSF" (merged group)
"""

import argparse
import csv
import math
import itertools
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--diversity",     nargs="+", required=True,
                        help="Per-sample *_srna_diversity.tsv files")
    parser.add_argument("--metrics",       nargs="+", required=True,
                        help="Per-sample *_diversity_metrics.tsv files")
    parser.add_argument("--out_diversity", default="all_samples_srna_diversity.tsv")
    parser.add_argument("--out_alpha",     default="all_samples_alpha_diversity.tsv")
    parser.add_argument("--out_beta",      default="all_samples_beta_diversity.tsv")
    parser.add_argument("--pseudocount",   type=float, default=1e-6,
                        help="Pseudocount for Aitchison CLR transform")
    parser.add_argument("--w_seg",         default=None,
                        help="Comma-separated segment tokens for per-group Whittaker "
                             "beta. Pipe-separated aliases within a token are merged "
                             "into one group. E.g. 'AWF,CL,LSW|LSF'")
    return parser.parse_args()


def parse_w_seg(w_seg_str):
    """
    Parse the --w_seg argument into a list of (group_label, [alias, ...]) tuples.

    Input:  "AWF,CL,LSW|LSF"
    Output: [
                ("AWF",     ["AWF"]),
                ("CL",      ["CL"]),
                ("LSW|LSF", ["LSW", "LSF"]),
            ]

    The group_label is used as the metric name in the output TSV.
    A sample matches a group if its sample_id contains ANY of the aliases
    (case-sensitive).
    """
    if not w_seg_str:
        return []

    groups = []
    for token in w_seg_str.split(","):
        token   = token.strip()
        aliases = [a.strip() for a in token.split("|") if a.strip()]
        if aliases:
            label = "|".join(aliases)
            groups.append((label, aliases))
    return groups


def sample_matches_group(sample_id, aliases):
    """Return True if sample_id contains any of the alias strings."""
    return any(alias in sample_id for alias in aliases)


# -----------------------------------------------------------------------------
# File loading
# -----------------------------------------------------------------------------

def merge_tsv_files(file_list, out_path):
    """Concatenate per-sample TSV files with a shared header into one file."""
    header_written = False
    with open(out_path, "w", newline="\n") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        for tsv_path in sorted(file_list):
            with open(tsv_path, newline="") as in_f:
                reader = csv.reader(in_f, delimiter="\t")
                header = next(reader)
                if not header_written:
                    writer.writerow(header)
                    header_written = True
                for row in reader:
                    writer.writerow(row)
    print(f"Written: {out_path}")


def load_proportions(diversity_files):
    """
    Load per-sample diversity TSVs.
    Returns:
        proportions  : { sample_id: { feature_type: fraction } }
        all_features : sorted list of all feature types across all samples
    """
    proportions  = defaultdict(dict)
    all_features = set()

    for tsv_path in diversity_files:
        with open(tsv_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample  = row["sample_id"].strip()
                feature = row["feature_type"].strip()
                frac    = float(row["fraction"])
                proportions[sample][feature] = frac
                all_features.add(feature)

    all_features = sorted(all_features)
    for sample in proportions:
        for feature in all_features:
            proportions[sample].setdefault(feature, 0.0)

    return dict(proportions), all_features


def load_alpha_metrics(metrics_files):
    """
    Load per-sample metrics TSVs.
    Returns: { sample_id: { metric: value } }
    """
    alpha = defaultdict(dict)
    for tsv_path in metrics_files:
        with open(tsv_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row["sample_id"] == "sample_id":
                    continue
                alpha[row["sample_id"]][row["metric"]] = row["value"]
    return dict(alpha)


# -----------------------------------------------------------------------------
# Beta diversity metrics
# -----------------------------------------------------------------------------

def bray_curtis(p1, p2, features):
    """
    Bray-Curtis dissimilarity.
    BC = sum(|p1k - p2k|) / sum(p1k + p2k)   [0–1]
    """
    numerator   = sum(abs(p1[f] - p2[f]) for f in features)
    denominator = sum(p1[f] + p2[f]      for f in features)
    return numerator / denominator if denominator > 0 else 0.0


def clr_transform(proportions, features, pseudocount):
    """Centered log-ratio transformation for compositional data."""
    log_vals = [math.log(proportions[f] + pseudocount) for f in features]
    mean_log = sum(log_vals) / len(log_vals)
    return {f: log_vals[i] - mean_log for i, f in enumerate(features)}


def aitchison_distance(p1, p2, features, pseudocount):
    """
    Aitchison distance: Euclidean distance in CLR-transformed space.
    Statistically rigorous for compositional data.   [0 → unbounded]
    """
    clr1 = clr_transform(p1, features, pseudocount)
    clr2 = clr_transform(p2, features, pseudocount)
    return math.sqrt(sum((clr1[f] - clr2[f]) ** 2 for f in features))


def whittaker_beta(proportions, all_features, sample_subset=None):
    """
    Whittaker's beta diversity.
    beta_W = gamma / mean(alpha) - 1

    gamma       = feature types observed in at least one sample in the subset
    mean(alpha) = mean number of feature types per sample in the subset

    If sample_subset is None, all samples are used.
    Returns (beta_W, n_samples_used) so the caller can record group size.
    """
    samples = sample_subset if sample_subset is not None else list(proportions.keys())

    if not samples:
        return 0.0, 0

    gamma      = sum(
        1 for f in all_features
        if any(proportions[s][f] > 0 for s in samples)
    )
    alpha_vals = [
        sum(1 for f in all_features if proportions[s][f] > 0)
        for s in samples
    ]
    alpha_mean = sum(alpha_vals) / len(alpha_vals) if alpha_vals else 0
    beta       = (gamma / alpha_mean) - 1 if alpha_mean > 0 else 0.0

    return beta, len(samples)


# -----------------------------------------------------------------------------
# Output writers
# -----------------------------------------------------------------------------

def write_alpha(alpha_metrics, samples, out_path):
    """Write wide-format alpha diversity table — one row per sample."""
    metric_order = [
        "total_library_reads",
        "total_mapped_reads",
        "mapping_rate",
        "feature_types_detected",
        "feature_types",
        "shannon_entropy_H",
        "max_possible_H_log2n",
        "normalized_H",
        "simpson_diversity_1_D",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id"] + metric_order)
        for sample in sorted(samples):
            row = [sample] + [
                alpha_metrics.get(sample, {}).get(m, "NA")
                for m in metric_order
            ]
            writer.writerow(row)
    print(f"Written: {out_path}")


def write_beta(proportions, all_features, pseudocount, w_seg_groups, out_path):
    """
    Write beta diversity output with three sections:
      1. Whittaker summary (global + per w_seg group if provided)
      2. Bray-Curtis pairwise matrix (all samples)
      3. Aitchison pairwise distance matrix (all samples)
    """
    samples   = sorted(proportions.keys())

    # --- Whittaker: global ---
    beta_global, n_global = whittaker_beta(proportions, all_features)

    # --- Whittaker: per segment group ---
    seg_results = []
    for label, aliases in w_seg_groups:
        subset = [s for s in samples if sample_matches_group(s, aliases)]
        beta_seg, n_seg = whittaker_beta(proportions, all_features,
                                         sample_subset=subset)
        seg_results.append((label, aliases, beta_seg, n_seg))
        print(f"Whittaker [{label}]: {beta_seg:.4f}  (n={n_seg})")

    # --- Pairwise matrices ---
    bc_matrix  = {s: {} for s in samples}
    ait_matrix = {s: {} for s in samples}

    for s1, s2 in itertools.product(samples, samples):
        bc_matrix[s1][s2]  = bray_curtis(
            proportions[s1], proportions[s2], all_features
        )
        ait_matrix[s1][s2] = aitchison_distance(
            proportions[s1], proportions[s2],
            all_features, pseudocount
        )

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Section 1: Whittaker summary
        writer.writerow(["## Beta Diversity Summary"])
        writer.writerow(["metric", "value", "n_samples"])
        writer.writerow(["n_samples_global",          str(len(samples)),  ""])
        writer.writerow(["n_feature_types_global",    str(len(all_features)), ""])
        writer.writerow(["feature_types_global",      ",".join(all_features), ""])
        writer.writerow(["whittaker_beta_global",     f"{beta_global:.4f}",
                          str(n_global)])

        # Per-segment Whittaker rows (only present if --w_seg was supplied)
        for label, aliases, beta_seg, n_seg in seg_results:
            writer.writerow([
                f"whittaker_beta_{label}",
                f"{beta_seg:.4f}",
                str(n_seg)
            ])

        writer.writerow([])

        # Section 2: Bray-Curtis matrix
        writer.writerow(["## Bray-Curtis Dissimilarity Matrix"])
        writer.writerow(["sample_id"] + samples)
        for s1 in samples:
            writer.writerow(
                [s1] + [f"{bc_matrix[s1][s2]:.6f}" for s2 in samples]
            )
        writer.writerow([])

        # Section 3: Aitchison distance matrix
        writer.writerow(["## Aitchison Distance Matrix"])
        writer.writerow(["sample_id"] + samples)
        for s1 in samples:
            writer.writerow(
                [s1] + [f"{ait_matrix[s1][s2]:.6f}" for s2 in samples]
            )

    print(f"Written: {out_path}")
    print(f"Whittaker beta (global): {beta_global:.4f}  (n={n_global})")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args         = parse_args()
    w_seg_groups = parse_w_seg(args.w_seg)

    if w_seg_groups:
        labels = [label for label, _ in w_seg_groups]
        print(f"Segment groups for Whittaker beta: {labels}")
    else:
        print("No --w_seg provided. Global Whittaker beta only.")

    # Merge per-feature count tables unchanged
    merge_tsv_files(args.diversity, args.out_diversity)

    # Load data
    proportions, all_features = load_proportions(args.diversity)
    alpha_metrics             = load_alpha_metrics(args.metrics)
    samples                   = sorted(proportions.keys())

    # Write outputs
    write_alpha(alpha_metrics, samples, args.out_alpha)
    write_beta(proportions, all_features, args.pseudocount,
               w_seg_groups, args.out_beta)


if __name__ == "__main__":
    main()