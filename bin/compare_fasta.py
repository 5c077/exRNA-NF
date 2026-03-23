#!/usr/bin/env python3
"""
compare_fasta.py

Compares two FASTA files to assess sequence similarity. Designed as a
sanity check to verify that pipeline-generated tRNA/rRNA FASTAs are
comparable to public database releases.

Usage:
    python3 compare_fasta.py --query pipeline.fa --reference database.fa
    python3 compare_fasta.py --query Ath_tRNA.fa --reference GtRNAdb_Ath.fa --out report.tsv
    python3 compare_fasta.py --query Ath_rRNA.fa --reference SILVA_Ath.fa --identity 0.90

Metrics reported:
    - Sequence count comparison
    - Length distribution comparison
    - Exact sequence matches
    - Identity-based matches (via pairwise alignment)
    - k-mer based similarity (fast, no alignment required)
    - Sequence type composition (for tRNA: anticodon/amino acid breakdown)
"""

import argparse
import sys
import math
from collections import defaultdict, Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare two FASTA files for sanity checking pipeline outputs"
    )
    parser.add_argument("--query",      required=True,
                        help="Pipeline-generated FASTA (tRNA or rRNA)")
    parser.add_argument("--reference",  required=True,
                        help="Public database FASTA to compare against")
    parser.add_argument("--out",        default=None,
                        help="Output TSV report path (default: stdout)")
    parser.add_argument("--identity",   type=float, default=0.95,
                        help="Minimum sequence identity for a match [default: 0.95]")
    parser.add_argument("--kmer_size",  type=int, default=8,
                        help="k-mer size for sequence similarity [default: 8]")
    parser.add_argument("--align",      action="store_true",
                        help="Perform pairwise alignment for identity scoring "
                             "(slow for large FASTAs, skip for >10k sequences)")
    parser.add_argument("--max_align",  type=int, default=500,
                        help="Max sequences to align per file if --align is set "
                             "[default: 500]")
    return parser.parse_args()


# -----------------------------------------------------------------------------
# FASTA parsing
# -----------------------------------------------------------------------------

def parse_fasta(path):
    """
    Parse a FASTA file into an ordered list of (header, sequence) tuples.
    Returns: list of (header, seq), dict of {seq: [headers]} for exact matching.
    """
    records    = []
    seq_index  = defaultdict(list)   # seq -> list of headers (handles duplicates)
    header     = None
    seq_parts  = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).upper()
                    records.append((header, seq))
                    seq_index[seq].append(header)
                header    = line[1:].split()[0]   # truncate at first whitespace
                seq_parts = []
            elif line:
                seq_parts.append(line)

    if header is not None:
        seq = "".join(seq_parts).upper()
        records.append((header, seq))
        seq_index[seq].append(header)

    return records, dict(seq_index)


# -----------------------------------------------------------------------------
# Length distribution
# -----------------------------------------------------------------------------

def length_stats(records):
    """Compute summary statistics for sequence lengths."""
    if not records:
        return {}
    lengths = sorted(len(seq) for _, seq in records)
    n       = len(lengths)
    total   = sum(lengths)
    mean    = total / n
    median  = lengths[n // 2] if n % 2 == 1 else \
              (lengths[n // 2 - 1] + lengths[n // 2]) / 2
    variance = sum((l - mean) ** 2 for l in lengths) / n
    sd       = math.sqrt(variance)

    return {
        "count":  n,
        "min":    lengths[0],
        "max":    lengths[-1],
        "mean":   round(mean, 2),
        "median": median,
        "sd":     round(sd, 2),
        "total_nt": total,
    }


def length_histogram(records, bins=10):
    """Return a simple text histogram of sequence lengths."""
    if not records:
        return []
    lengths  = [len(seq) for _, seq in records]
    min_l    = min(lengths)
    max_l    = max(lengths)
    if min_l == max_l:
        return [(f"{min_l}", len(lengths))]
    bin_size = max(1, (max_l - min_l) // bins)
    hist     = defaultdict(int)
    for l in lengths:
        bucket = (l - min_l) // bin_size
        hist[bucket] += 1
    rows = []
    for b in sorted(hist):
        lo = min_l + b * bin_size
        hi = lo + bin_size - 1
        rows.append((f"{lo}–{hi}", hist[b]))
    return rows


# -----------------------------------------------------------------------------
# Exact match analysis
# -----------------------------------------------------------------------------

def exact_match_analysis(query_records, query_index,
                          ref_records,   ref_index):
    """
    Find sequences present in both files as exact matches.
    Handles the case where sequences may appear in one but not the other.
    """
    query_seqs = set(query_index.keys())
    ref_seqs   = set(ref_index.keys())

    exact_both     = query_seqs & ref_seqs
    query_only     = query_seqs - ref_seqs
    ref_only       = ref_seqs   - query_seqs

    return {
        "exact_match_count":        len(exact_both),
        "query_only_count":         len(query_only),
        "ref_only_count":           len(ref_only),
        "pct_query_in_ref":         round(len(exact_both) / len(query_seqs) * 100, 2)
                                    if query_seqs else 0.0,
        "pct_ref_in_query":         round(len(exact_both) / len(ref_seqs) * 100, 2)
                                    if ref_seqs else 0.0,
    }


# -----------------------------------------------------------------------------
# k-mer similarity
# -----------------------------------------------------------------------------

def build_kmer_profile(seq, k):
    """Return Counter of all k-mers in a sequence."""
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))


def jaccard_kmer(profile1, profile2):
    """
    Jaccard similarity between two k-mer profiles.
    J = |intersection| / |union|   [0–1]
    """
    keys1 = set(profile1.keys())
    keys2 = set(profile2.keys())
    inter = sum(min(profile1[k], profile2[k]) for k in keys1 & keys2)
    union = sum(profile1.values()) + sum(profile2.values()) - inter
    return inter / union if union > 0 else 0.0


def kmer_similarity_analysis(query_records, ref_records, k):
    """
    Compute k-mer based similarity between the two FASTA files.

    Builds a combined k-mer profile for all sequences in each file
    and computes Jaccard similarity. This is fast and alignment-free,
    useful for large rRNA files.
    """
    print(f"  Building {k}-mer profiles...", file=sys.stderr)

    q_profile = Counter()
    r_profile = Counter()

    for _, seq in query_records:
        q_profile.update(build_kmer_profile(seq, k))

    for _, seq in ref_records:
        r_profile.update(build_kmer_profile(seq, k))

    jaccard = jaccard_kmer(q_profile, r_profile)

    # Per-sequence best-match k-mer similarity
    # Sample up to 200 query sequences to keep runtime reasonable
    sample_q = query_records[:200]
    sample_r = ref_records[:200]

    ref_profiles = [(h, build_kmer_profile(seq, k)) for h, seq in sample_r]

    per_seq_scores = []
    for q_hdr, q_seq in sample_q:
        q_prof = build_kmer_profile(q_seq, k)
        if not ref_profiles:
            continue
        best = max(jaccard_kmer(q_prof, r_prof)
                   for _, r_prof in ref_profiles)
        per_seq_scores.append(best)

    mean_best = round(sum(per_seq_scores) / len(per_seq_scores), 4) \
                if per_seq_scores else 0.0

    return {
        "kmer_size":                k,
        "jaccard_global":           round(jaccard, 4),
        "mean_best_kmer_per_seq":   mean_best,
        "n_query_sampled":          len(sample_q),
    }


# -----------------------------------------------------------------------------
# Pairwise alignment identity (optional)
# -----------------------------------------------------------------------------

def needleman_wunsch_identity(seq1, seq2,
                               match=1, mismatch=-1, gap=-2):
    """
    Simple Needleman-Wunsch global alignment.
    Returns sequence identity as fraction of aligned positions.
    Only practical for sequences up to ~500 nt.
    """
    n, m = len(seq1), len(seq2)

    # Limit to prevent memory issues on long sequences
    if n > 600 or m > 600:
        # Fall back to k-mer proxy for long sequences
        k = 8
        p1 = build_kmer_profile(seq1, k)
        p2 = build_kmer_profile(seq2, k)
        return jaccard_kmer(p1, p2)

    # Score matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = gap * i
    for j in range(m + 1):
        dp[0][j] = gap * j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            dp[i][j] = max(
                dp[i - 1][j - 1] + s,
                dp[i - 1][j]     + gap,
                dp[i][j - 1]     + gap,
            )

    # Traceback to count matches
    i, j     = n, m
    matches  = 0
    aligned  = 0

    while i > 0 and j > 0:
        s = match if seq1[i - 1] == seq2[j - 1] else mismatch
        if dp[i][j] == dp[i - 1][j - 1] + s:
            if seq1[i - 1] == seq2[j - 1]:
                matches += 1
            i -= 1
            j -= 1
            aligned += 1
        elif dp[i][j] == dp[i - 1][j] + gap:
            i -= 1
            aligned += 1
        else:
            j -= 1
            aligned += 1

    return matches / aligned if aligned > 0 else 0.0


def alignment_identity_analysis(query_records, ref_records,
                                 min_identity, max_seqs):
    """
    For each query sequence, find the best-identity match in the reference
    using pairwise global alignment.

    Limited to max_seqs sequences to keep runtime reasonable.
    Reports the distribution of best-match identity scores.
    """
    sample_q = query_records[:max_seqs]
    sample_r = ref_records[:max_seqs]

    print(f"  Aligning {len(sample_q)} query vs {len(sample_r)} reference "
          f"sequences...", file=sys.stderr)

    above_threshold = 0
    identity_scores = []

    for i, (q_hdr, q_seq) in enumerate(sample_q):
        if i % 50 == 0:
            print(f"    {i}/{len(sample_q)} aligned", file=sys.stderr)

        best = 0.0
        for _, r_seq in sample_r:
            identity = needleman_wunsch_identity(q_seq, r_seq)
            if identity > best:
                best = identity
            if best == 1.0:
                break   # exact match found, stop searching

        identity_scores.append(best)
        if best >= min_identity:
            above_threshold += 1

    mean_id = round(sum(identity_scores) / len(identity_scores), 4) \
              if identity_scores else 0.0

    bins = [0.0, 0.7, 0.8, 0.9, 0.95, 0.99, 1.01]
    labels = ["<0.70", "0.70–0.80", "0.80–0.90",
              "0.90–0.95", "0.95–0.99", "1.00"]
    dist = Counter()
    for score in identity_scores:
        for b in range(len(bins) - 1):
            if bins[b] <= score < bins[b + 1]:
                dist[labels[b]] += 1
                break

    return {
        "n_query_aligned":          len(sample_q),
        "n_ref_aligned":            len(sample_r),
        "mean_best_identity":       mean_id,
        "pct_above_threshold":      round(above_threshold /
                                         len(sample_q) * 100, 2)
                                    if sample_q else 0.0,
        "identity_threshold":       min_identity,
        "identity_distribution":    dist,
    }


# -----------------------------------------------------------------------------
# tRNA-specific analysis
# -----------------------------------------------------------------------------

def parse_trna_anticodon(header):
    """
    Attempt to extract amino acid and anticodon from tRNA headers.
    Handles common formats:
      tRNAscan-SE: Chr1.trna1-AlaAGC  ->  Ala, AGC
      GtRNAdb:     Arabidopsis_thaliana_chr1.trna1-AlaAGC
      Generic:     looks for 3-letter AA code followed by 3-letter anticodon
    """
    import re
    # Pattern: 3-letter amino acid + 3-letter anticodon
    match = re.search(r'([A-Z][a-z]{2})([ACGT]{3})', header)
    if match:
        return match.group(1), match.group(2)
    # Suppressor tRNAs
    if "Sup" in header or "sup" in header:
        return "Sup", "UNK"
    return "Unknown", "Unknown"


def trna_composition(records):
    """
    Break down tRNA records by amino acid and anticodon.
    Returns Counter of amino acids and Counter of anticodons.
    """
    aa_counts       = Counter()
    anticodon_counts = Counter()

    for header, seq in records:
        aa, anticodon = parse_trna_anticodon(header)
        aa_counts[aa]              += 1
        anticodon_counts[anticodon] += 1

    return aa_counts, anticodon_counts


# -----------------------------------------------------------------------------
# Report writing
# -----------------------------------------------------------------------------

def write_report(args, query_records, ref_records,
                 exact_results, kmer_results,
                 align_results, out):
    """Write the full comparison report as a structured TSV."""

    q_stats = length_stats(query_records)
    r_stats = length_stats(ref_records)

    lines = []

    def row(*cols):
        lines.append("\t".join(str(c) for c in cols))

    # Header
    row("## FASTA Comparison Report")
    row("query",     args.query)
    row("reference", args.reference)
    row("")

    # Section 1: Sequence counts
    row("## Sequence Counts")
    row("metric",                   "query",            "reference")
    row("sequence_count",           q_stats["count"],   r_stats["count"])
    row("total_nucleotides",        q_stats["total_nt"],r_stats["total_nt"])
    row("")

    # Section 2: Length distribution
    row("## Length Distribution")
    row("metric",   "query",            "reference")
    row("min_nt",   q_stats["min"],     r_stats["min"])
    row("max_nt",   q_stats["max"],     r_stats["max"])
    row("mean_nt",  q_stats["mean"],    r_stats["mean"])
    row("median_nt",q_stats["median"],  r_stats["median"])
    row("sd_nt",    q_stats["sd"],      r_stats["sd"])
    row("")

    # Section 3: Exact match analysis
    row("## Exact Sequence Matching")
    row("metric",                       "value")
    row("exact_matches",                exact_results["exact_match_count"])
    row("query_only_sequences",         exact_results["query_only_count"])
    row("reference_only_sequences",     exact_results["ref_only_count"])
    row("pct_query_found_in_ref",       f"{exact_results['pct_query_in_ref']}%")
    row("pct_ref_found_in_query",       f"{exact_results['pct_ref_in_query']}%")
    row("")

    # Section 4: k-mer similarity
    row("## k-mer Similarity")
    row("metric",                       "value")
    row("kmer_size",                    kmer_results["kmer_size"])
    row("jaccard_global",               kmer_results["jaccard_global"])
    row("mean_best_kmer_per_seq",       kmer_results["mean_best_kmer_per_seq"])
    row("n_query_sequences_sampled",    kmer_results["n_query_sampled"])
    row("")

    # Section 5: Alignment identity (if run)
    if align_results:
        row("## Pairwise Alignment Identity")
        row("metric",                   "value")
        row("n_query_aligned",          align_results["n_query_aligned"])
        row("n_ref_aligned",            align_results["n_ref_aligned"])
        row("mean_best_identity",       align_results["mean_best_identity"])
        row("pct_above_threshold",
            f"{align_results['pct_above_threshold']}% "
            f"(threshold={align_results['identity_threshold']})")
        row("")
        row("## Identity Distribution")
        row("identity_bin",             "count")
        for label, count in sorted(align_results["identity_distribution"].items()):
            row(label, count)
        row("")

    # Section 6: tRNA-specific composition (if headers suggest tRNA content)
    sample_headers = " ".join(h for h, _ in query_records[:20])
    is_trna = any(kw in sample_headers.lower()
                  for kw in ["trna", "ala", "gly", "val", "leu",
                              "ile", "pro", "phe", "trp", "met"])
    if is_trna:
        row("## tRNA Composition (query)")
        q_aa, q_ac = trna_composition(query_records)
        r_aa, r_ac = trna_composition(ref_records)

        all_aa = sorted(set(q_aa.keys()) | set(r_aa.keys()))
        row("amino_acid",   "query_count",  "ref_count",    "difference")
        for aa in all_aa:
            diff = q_aa[aa] - r_aa[aa]
            row(aa, q_aa[aa], r_aa[aa], diff)
        row("")

        row("## tRNA Anticodon Comparison")
        all_ac = sorted(set(q_ac.keys()) | set(r_ac.keys()))
        row("anticodon",    "query_count",  "ref_count",    "difference")
        for ac in all_ac:
            diff = q_ac[ac] - r_ac[ac]
            row(ac, q_ac[ac], r_ac[ac], diff)
        row("")

    # Section 7: Interpretation summary
    row("## Interpretation")
    jaccard     = kmer_results["jaccard_global"]
    pct_q_in_r  = exact_results["pct_query_in_ref"]

    if pct_q_in_r >= 90:
        exact_verdict = "PASS — >90% of query sequences found exactly in reference"
    elif pct_q_in_r >= 70:
        exact_verdict = "WARN — 70–90% exact match; minor version differences expected"
    else:
        exact_verdict = "FAIL — <70% exact match; significant divergence detected"

    if jaccard >= 0.80:
        kmer_verdict = "PASS — high k-mer similarity (Jaccard ≥ 0.80)"
    elif jaccard >= 0.60:
        kmer_verdict = "WARN — moderate k-mer similarity (0.60–0.80)"
    else:
        kmer_verdict = "FAIL — low k-mer similarity (< 0.60); check organism/version"

    row("exact_match_verdict",  exact_verdict)
    row("kmer_verdict",         kmer_verdict)

    output_str = "\n".join(lines)

    if out:
        with open(out, "w") as f:
            f.write(output_str + "\n")
        print(f"Report written: {out}", file=sys.stderr)
    else:
        print(output_str)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()

    print(f"Loading query:     {args.query}",     file=sys.stderr)
    query_records, query_index = parse_fasta(args.query)
    print(f"  {len(query_records)} sequences",    file=sys.stderr)

    print(f"Loading reference: {args.reference}", file=sys.stderr)
    ref_records, ref_index = parse_fasta(args.reference)
    print(f"  {len(ref_records)} sequences",      file=sys.stderr)

    if not query_records:
        print("ERROR: query FASTA is empty", file=sys.stderr)
        sys.exit(1)
    if not ref_records:
        print("ERROR: reference FASTA is empty", file=sys.stderr)
        sys.exit(1)

    # Exact match analysis
    print("Running exact match analysis...", file=sys.stderr)
    exact_results = exact_match_analysis(
        query_records, query_index, ref_records, ref_index
    )

    # k-mer similarity
    print(f"Running {args.kmer_size}-mer similarity analysis...", file=sys.stderr)
    kmer_results = kmer_similarity_analysis(
        query_records, ref_records, args.kmer_size
    )

    # Pairwise alignment (optional)
    align_results = None
    if args.align:
        print("Running pairwise alignment identity analysis...", file=sys.stderr)
        align_results = alignment_identity_analysis(
            query_records, ref_records,
            args.identity, args.max_align
        )

    # Write report
    write_report(args, query_records, ref_records,
                 exact_results, kmer_results,
                 align_results, args.out)


if __name__ == "__main__":
    main()