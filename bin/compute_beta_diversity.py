#!/usr/bin/env python3
"""
Compute beta diversity metrics across all sRNA library samples.
Input: all_samples_srna_diversity.tsv produced by MERGE_DIVERSITY_OUTPUTS
Output:
  - bray_curtis_matrix.tsv    pairwise Bray-Curtis dissimilarity matrix
  - aitchison_matrix.tsv      pairwise Aitchison distance matrix
  - beta_diversity_summary.tsv whittaker beta + per-sample summaries
"""

import argparse
import csv
import math
import itertools
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--diversity", required=True,
                        help="all_samples_srna_diversity.tsv")
    parser.add_argument("--out_bray_curtis",  default="bray_curtis_matrix.tsv")
    parser.add_argument("--out_aitchison",    default="aitchison_matrix.tsv")
    parser.add_argument("--out_summary",      default="beta_diversity_summary.tsv")
    parser.add_argument("--pseudocount",      type=float, default=1e-6,
                        help="Pseudocount added before log transform to handle zeros")
    return parser.parse_args()


# -----------------------------------------------------------------------------
# Load diversity table
# -----------------------------------------------------------------------------

def load_proportions(diversity_path):
    """
    Load all_samples_srna_diversity.tsv and return a dict of:
        { sample_id: { feature_type: fraction } }
    Uses the 'fraction' column which already sums to 1 per sample.
    """
    proportions = defaultdict(dict)
    all_features = set()

    with open(diversity_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample  = row["sample_id"]
            feature = row["feature_type"]
            frac    = float(row["fraction"])
            proportions[sample][feature] = frac
            all_features.add(feature)

    return dict(proportions), sorted(all_features)


def fill_missing_features(proportions, all_features):
    """
    Ensure every sample has an entry for every feature type.
    Missing features are filled with 0.
    """
    for sample in proportions:
        for feature in all_features:
            if feature not in proportions[sample]:
                proportions[sample][feature] = 0.0
    return proportions


# -----------------------------------------------------------------------------
# Beta diversity metrics
# -----------------------------------------------------------------------------

def bray_curtis(p1, p2, features):
    """
    Bray-Curtis dissimilarity between two samples.
    BC = sum(|p1k - p2k|) / sum(p1k + p2k)
    Range: 0 (identical) to 1 (completely different)
    """
    numerator   = sum(abs(p1[f] - p2[f]) for f in features)
    denominator = sum(p1[f] + p2[f] for f in features)
    return numerator / denominator if denominator > 0 else 0.0


def clr_transform(proportions, features, pseudocount):
    """
    Centered log-ratio transformation for compositional data.
    clr(p_k) = log(p_k + pseudo) - mean(log(p_k + pseudo))
    """
    log_values = [math.log(proportions[f] + pseudocount) for f in features]
    mean_log   = sum(log_values) / len(log_values)
    return {f: log_values[i] - mean_log for i, f in enumerate(features)}


def aitchison_distance(p1, p2, features, pseudocount):
    """
    Aitchison distance: Euclidean distance in clr-transformed space.
    More appropriate than Bray-Curtis for compositional data because
    it respects the simplex geometry of proportions.
    Range: 0 (identical composition) to unbounded positive values.
    """
    clr1 = clr_transform(p1, features, pseudocount)
    clr2 = clr_transform(p2, features, pseudocount)
    return math.sqrt(sum((clr1[f] - clr2[f]) ** 2 for f in features))


def whittaker_beta(proportions, all_features):
    """
    Whittaker's beta diversity: beta_W = gamma/mean(alpha) - 1
    gamma = total feature types observed across all samples
    alpha = feature types observed per sample (fraction > 0)
    Range: 0 (all samples identical) to n_samples - 1 (no shared features)
    """
    gamma      = sum(1 for f in all_features
                     if any(proportions[s][f] > 0 for s in proportions))
    alpha_vals = [sum(1 for f in all_features if proportions[s][f] > 0)
                  for s in proportions]
    alpha_mean = sum(alpha_vals) / len(alpha_vals) if alpha_vals else 0
    return (gamma / alpha_mean) - 1 if alpha_mean > 0 else 0.0


# -----------------------------------------------------------------------------
# Write output matrices
# -----------------------------------------------------------------------------

def write_distance_matrix(matrix, samples, out_path):
    """Write a symmetric pairwise distance matrix as TSV."""
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id"] + samples)
        for s1 in samples:
            row = [s1] + [f"{matrix[s1][s2]:.6f}" for s2 in samples]
            writer.writerow(row)


def write_summary(proportions, all_features, whittaker, out_path):
    """Write beta diversity summary including Whittaker and per-sample feature counts."""
    samples = sorted(proportions.keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["n_samples",              str(len(samples))])
        writer.writerow(["n_feature_types_global", str(len(all_features))])
        writer.writerow(["feature_types_global",   ",".join(all_features)])
        writer.writerow(["whittaker_beta",          f"{whittaker:.4f}"])
        writer.writerow([])
        writer.writerow(["sample_id", "n_feature_types_detected",
                         "feature_types_detected"])
        for s in samples:
            detected = [f for f in all_features if proportions[s][f] > 0]
            writer.writerow([s, len(detected), ",".join(detected)])


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()

    proportions, all_features = load_proportions(args.diversity)
    proportions = fill_missing_features(proportions, all_features)
    samples     = sorted(proportions.keys())

    print(f"Loaded {len(samples)} samples x {len(all_features)} feature types")

    # Compute pairwise Bray-Curtis matrix
    bc_matrix = {s: {} for s in samples}
    for s1, s2 in itertools.product(samples, samples):
        bc_matrix[s1][s2] = bray_curtis(
            proportions[s1], proportions[s2], all_features
        )

    # Compute pairwise Aitchison distance matrix
    ait_matrix = {s: {} for s in samples}
    for s1, s2 in itertools.product(samples, samples):
        ait_matrix[s1][s2] = aitchison_distance(
            proportions[s1], proportions[s2],
            all_features, args.pseudocount
        )

    # Compute Whittaker beta
    whittaker = whittaker_beta(proportions, all_features)

    # Write outputs
    write_distance_matrix(bc_matrix,  samples, args.out_bray_curtis)
    write_distance_matrix(ait_matrix, samples, args.out_aitchison)
    write_summary(proportions, all_features, whittaker, args.out_summary)

    print(f"Whittaker beta diversity: {whittaker:.4f}")
    print(f"Bray-Curtis matrix  -> {args.out_bray_curtis}")
    print(f"Aitchison matrix    -> {args.out_aitchison}")
    print(f"Summary             -> {args.out_summary}")


if __name__ == "__main__":
    main()