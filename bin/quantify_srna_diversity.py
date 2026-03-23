#!/usr/bin/env python3

import pysam
import argparse, sys
import csv
import math
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id",   required=True)
    parser.add_argument("--bam",         required=True)
    parser.add_argument("--combined_fa", required=True,
                        help="Labeled combined FASTA used to build the bowtie index. "
                             "Headers must follow the format: >featureType|sequenceId")
    parser.add_argument("--bowtie2_stats", required=True)
    parser.add_argument("--out",         default="srna_diversity.tsv")
    parser.add_argument("--metrics",     default="diversity_metrics.tsv")
    return parser.parse_args()

# -----------------------------------------------------------------------------
# Get mapping statistics from the Bowtie2 --stdout--
# -----------------------------------------------------------------------------
def parse_bowtie2_total_reads(stats_path):
    """
    Parse total input reads from bowtie2 stderr stats file.
    bowtie2 reports: 'X reads; of these:'
    """
    try:
        with open(stats_path) as f:
            content = f.read()
            # Print first few lines for debugging
            print(f"[DEBUG] bowtie2 stats content:\n{content[:200]}",
                  file=sys.stderr)
            for line in content.splitlines():
                line = line.strip()
                if "reads; of these:" in line:
                    total = int(line.split()[0])
                    print(f"[DEBUG] Parsed total reads: {total}",
                          file=sys.stderr)
                    return total
        print(f"[WARN] Could not find 'reads; of these:' in {stats_path}",
              file=sys.stderr)
    except FileNotFoundError:
        print(f"[WARN] bowtie2 stats file not found: {stats_path}",
              file=sys.stderr)
    except Exception as e:
        print(f"[WARN] Could not parse bowtie2 stats {stats_path}: {e}",
              file=sys.stderr)
    return 0

# -----------------------------------------------------------------------------
# Label map from combined FASTA
# -----------------------------------------------------------------------------

def build_label_map(combined_fa_path):
    label_map = {}
    with open(combined_fa_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                # Truncate at first whitespace to match bowtie2 BAM reference naming
                header = header.split()[0]
                feature_type = header.split("|")[0]
                label_map[header] = feature_type
    print(f"Label map built: {len(label_map)} sequences across "
          f"{len(set(label_map.values()))} feature types")
    return label_map


# -----------------------------------------------------------------------------
# Read counting
# -----------------------------------------------------------------------------

def get_collapsed_count(read_name):
    """Recover abundance from fastx_collapser names (seq_N_xCOUNT)."""
    try:
        return int(read_name.split("_x")[-1])
    except (ValueError, IndexError):
        return 1


def collect_read_hits(bam_path, label_map):
    """
    Single-pass over the combined BAM.
    Records which feature types each read mapped to, and its abundance.
    A read mapping to miRNA|ath-miR166a and miRNA|ath-miR166b is counted
    as a single-feature read (both are miRNA), not a multimapper.
    """
    read_hits      = defaultdict(set)
    read_abundance = {}

    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = bam.references   # reference name list from BAM header

    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped:
            continue
        ref_name     = refs[aln.reference_id]
        feature_type = label_map.get(ref_name, "unknown")
        rname        = aln.query_name

        read_hits[rname].add(feature_type)
        if rname not in read_abundance:
            read_abundance[rname] = get_collapsed_count(rname)

    bam.close()
    return read_hits, read_abundance


# -----------------------------------------------------------------------------
# Fractional counting
# -----------------------------------------------------------------------------

def compute_fractional_counts(read_hits, read_abundance):
    feature_counts = defaultdict(float)
    feature_unique = defaultdict(float)
    total_reads    = 0.0

    for rname, hit_types in read_hits.items():
        abundance    = read_abundance[rname]
        total_reads += abundance
        share        = abundance / len(hit_types)

        for ft in hit_types:
            feature_counts[ft] += share
            if len(hit_types) == 1:
                feature_unique[ft] += abundance

    return feature_counts, feature_unique, total_reads


# -----------------------------------------------------------------------------
# Diversity metrics
# -----------------------------------------------------------------------------

def shannon_entropy(proportions):
    return -sum(p * math.log2(p) for p in proportions if p > 0)


def simpson_diversity(proportions):
    return 1.0 - sum(p ** 2 for p in proportions)


def rao_pd(proportions):
    """
    Rao's Quadratic Entropy with categorical distances (d_ij = 1 for i != j).
    Equivalent to Simpson's index under equal inter-class distances.
    Supply a distance matrix for true phylogenetic PD weighting.
    """
    return 1.0 - sum(p ** 2 for p in proportions)


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

def write_counts(sample_id, feature_counts, feature_unique, total_reads, total_library_reads, out_path):
    total_assigned = sum(feature_counts.values())
    rpm_scale      = 1_000_000 / total_library_reads if total_library_reads > 0 else 0.0

    with open(out_path, "w", newline="\n") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample_id", "feature_type",
            "fractional_count",
            "unique_feature_count",     # reads mapping to only this feature type
            "shared_feature_count",     # reads mapping to this + other feature types
            "RPM",
            "fraction",
            "percent"
        ])
        for ft, count in sorted(feature_counts.items(), key=lambda x: -x[1]):
            unique = feature_unique.get(ft, 0.0)
            shared = count - unique     # replaces multimapper_count
            frac   = count / total_assigned if total_assigned > 0 else 0.0
            writer.writerow([
                sample_id, ft,
                f"{count:.2f}",
                f"{unique:.2f}",
                f"{shared:.2f}",
                f"{count * rpm_scale:.3f}",
                f"{frac:.4f}",
                f"{frac * 100:.2f}"
            ])


def write_metrics(sample_id, feature_counts, total_reads, total_library_reads, metrics_path):
    total_assigned = sum(feature_counts.values())
    if total_assigned == 0:
        print(f"[WARN] No assigned reads for {sample_id}")
        return

    feature_list = sorted(feature_counts.keys())
    proportions  = [feature_counts[ft] / total_assigned for ft in feature_list]
    n_types      = len(feature_list)
    H            = shannon_entropy(proportions)
    H_max        = math.log2(n_types) if n_types > 1 else 0.0
    H_norm       = H / H_max if H_max > 0 else "NA"
    D            = simpson_diversity(proportions)

    with open(metrics_path, "w", newline="\n") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "metric", "value"])
        for metric, value in [
            ("total_library_reads",        f"{total_library_reads:.0f}"),
            ("total_mapped_reads",         f"{total_reads:.0f}"),
            ("mapping_rate",               f"{total_reads / total_library_reads:.4f}"
                                           if total_library_reads > 0 else "NA"),
            ("feature_types_detected",     str(n_types)),
            ("feature_types",              ",".join(feature_list)),
            ("shannon_entropy_H",          f"{H:.4f}"),
            ("max_possible_H_log2n",       f"{H_max:.4f}"),
            ("normalized_H",               f"{H_norm:.4f}"
                                           if isinstance(H_norm, float) else H_norm),
            ("simpson_diversity_1_D",      f"{D:.4f}"),
        ]:
            writer.writerow([sample_id, metric, value])

    print(f"[{sample_id}] H={H:.3f}  "
          f"H_norm={H_norm if isinstance(H_norm, str) else f'{H_norm:.3f}'}  "
          f"Simpson={D:.3f}  "
          f"Mapped={total_reads}  "
          f"Library={total_library_reads}  "
          f"Rate={total_reads/total_library_reads:.3f}"
          if total_library_reads > 0 else "")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args      = parse_args()
    label_map = build_label_map(args.combined_fa)

    # Parse total library size from bowtie2 stats
    total_library_reads = parse_bowtie2_total_reads(args.bowtie2_stats)

    read_hits, read_abundance = collect_read_hits(args.bam, label_map)
    feature_counts, feature_unique, total_reads = compute_fractional_counts(
        read_hits, read_abundance
    )

    write_counts(
        args.sample_id, feature_counts, feature_unique,
        total_reads, total_library_reads, args.out      # add total_library_reads
    )
    write_metrics(
        args.sample_id, feature_counts,
        total_reads, total_library_reads, args.metrics  # add total_library_reads
    )
    print(f"[{args.sample_id}] Done -> {args.out} | {args.metrics}")


if __name__ == "__main__":
    main()