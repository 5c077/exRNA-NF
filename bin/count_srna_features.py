#!/usr/bin/env python3
"""
count_srna_features.py

Per-feature-name fractional count and RPM quantification from the
combined annotation BAM. Produces one row per (library, feature_type,
feature_name) combination.

tRNA is split into three categories:
    tRNA_all   : all tRNA-mapping reads (10-50 nt)
    tRNA_tRF   : tRNA-derived fragments (10-27 nt, exclusive upper bound)
    tRNA_tiRNA : tRNA halves / tiRNAs   (28-50 nt, inclusive both bounds)

All other feature types are reported as-is with per-sequence-ID rows.

Usage:
    count_srna_features.py \
        --sample_id     <lib_name> \
        --bam           <combined_annotations.bam> \
        --combined_fa   <combined.fa> \
        --bowtie2_stats <bowtie2_stats.txt> \
        --out           <output.csv>
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# tRNA size class boundaries (exclusive upper / inclusive lower)
# tRF  : 10 <= length <= 27
# tiRNA: 28 <= length <= 50
# ---------------------------------------------------------------------------
TRF_MIN   = 10
TRF_MAX   = 27    # inclusive
TIRNA_MIN = 28    # inclusive
TIRNA_MAX = 50


def parse_args():
    parser = argparse.ArgumentParser(
        description="Per-feature-name fractional abundance counting"
    )
    parser.add_argument("--sample_id",     required=True)
    parser.add_argument("--bam",           required=True,
                        help="Combined annotation BAM")
    parser.add_argument("--genome_name",   required=True,        
                        help="Genome label derived from combined FASTA filename")
    parser.add_argument("--combined_fa",   required=True,
                        help="Labeled combined FASTA (featureType|seqId headers)")
    parser.add_argument("--bowtie2_stats", required=True,
                        help="bowtie2 stderr stats file for total library reads")
    parser.add_argument("--out",           required=True,
                        help="Output CSV file")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Library metadata parsing  (same logic as extractSizeDistribution)
# ---------------------------------------------------------------------------

def parse_library_metadata(sample_id, genome_name):
    """
    CHANGED: genome is now passed in directly from the combined_fa filename
    rather than parsed from sample_id, since lib_name from
    alignToCombinedAnnotations does not contain the _vs_genome_final suffix.

    Tissue is still parsed from the sample_id prefix.
    """
    # CHANGED: use genome_name parameter directly
    genome = genome_name

    # Tissue is still parseable from lib_name prefix
    tissue_tokens = ["AWF", "CL", "LSW", "LSF"]
    tissue        = "unknown"
    for token in tissue_tokens:
        if re.search(r'(?:^|_)' + token + r'(?:_|$)', sample_id):
            tissue = token
            break

    return genome, tissue


# ---------------------------------------------------------------------------
# Label map:  BAM ref_name  ->  (feature_type, feature_name)
# ---------------------------------------------------------------------------

def build_label_map_from_bam(bam_path):
    """
    CHANGED: build label map from BAM reference header instead of reading
    the combined FASTA.

    The combined FASTA for Lsa contains 8M+ sequences and causes a MemoryError
    when read into a Python dict. The BAM header only contains references that
    reads actually aligned to, which is orders of magnitude smaller.

    Reference names in the BAM follow the same featureType|sequenceId format
    as the FASTA headers, so no information is lost.
    """
    import pysam

    label_map = {}
    bam = pysam.AlignmentFile(bam_path, "rb")

    for ref_name in bam.references:
        # Truncate at whitespace (bowtie2 does this when writing the BAM)
        ref_name = ref_name.split()[0]
        parts    = ref_name.split("|", 1)
        if len(parts) == 2:
            feature_type, feature_name = parts
        else:
            feature_type = parts[0]
            feature_name = parts[0]
        label_map[ref_name] = (feature_type, feature_name)

    bam.close()

    n_types = len({ft for ft, _ in label_map.values()})
    print(f"Label map from BAM header: {len(label_map)} references across "
          f"{n_types} feature types", file=sys.stderr)
    return label_map

def enrich_label_map_with_full_headers(label_map, combined_fa_path):
    """
    Scan the combined FASTA and replace feature_name in label_map with
    the full original header line for sequences present in the BAM.

    Memory safe: only retains entries for sequences already in label_map,
    which is the subset of references that reads actually aligned to.
    Stops scanning as soon as all needed sequences are found.

    Combined FASTA header format:
        >featureType|sequenceId description tokens ...
    e.g.
        >miRNA|at1g01183 ath-MIR156a MI0000178 Arabidopsis thaliana miR156a
        >tRNA|Chr1.trna1-AlaAGC Ala (AGC) 3631-3820 Score: 88.4

    feature_name becomes everything after featureType| including description:
        at1g01183 ath-MIR156a MI0000178 Arabidopsis thaliana miR156a
        Chr1.trna1-AlaAGC Ala (AGC) 3631-3820 Score: 88.4
    """
    refs_needed = set(label_map.keys())
    enriched    = 0

    with open(combined_fa_path, errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                continue

            full_header = line[1:].strip()
            first_token = full_header.split()[0]    # what bowtie2 writes to BAM

            if first_token not in refs_needed:
                continue

            feature_type    = label_map[first_token][0]  # preserve type from BAM map
            rest_of_header  = full_header[len(first_token):].strip()

            # feature_name = sequenceId + everything after it in the header
            # e.g. "at1g01183 ath-MIR156a MI0000178 Arabidopsis thaliana miR156a"
            parts = first_token.split("|", 1)
            seq_id = parts[1] if len(parts) == 2 else first_token
            feature_name = f"{seq_id} {rest_of_header}".strip() if rest_of_header else seq_id

            label_map[first_token] = (feature_type, feature_name)
            enriched += 1

            # Stop scanning once all BAM references are found
            if enriched == len(refs_needed):
                break

    print(f"  Full headers enriched: {enriched}/{len(refs_needed)} references",
          file=sys.stderr)
    return label_map
# ---------------------------------------------------------------------------
# bowtie2 stats parsing
# ---------------------------------------------------------------------------

def parse_bowtie2_total_reads(stats_path):
    """Parse total input reads from bowtie2 stats file."""
    try:
        with open(stats_path) as f:
            for line in f:
                if "reads; of these:" in line:
                    return int(line.strip().split()[0])
    except Exception as e:
        print(f"[WARN] Could not parse bowtie2 stats {stats_path}: {e}",
              file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Read abundance (collapsed reads carry count in name)
# ---------------------------------------------------------------------------

def get_collapsed_count(read_name):
    """Recover abundance from fastx_collapser-style names (seq_N_xCOUNT)."""
    try:
        return int(read_name.split("_x")[-1])
    except (ValueError, IndexError):
        return 1


# ---------------------------------------------------------------------------
# Core counting logic
# ---------------------------------------------------------------------------

def collect_feature_counts(bam_path, label_map):
    """
    CHANGED: streaming name-sorted BAM approach to eliminate the
    read_hits dictionary that caused MemoryError on large libraries.

    Strategy:
        1. Sort BAM by query name into a temporary file using samtools
        2. Stream through name-sorted BAM one read group at a time
        3. For each group of alignments sharing a read name, compute
           the fractional share and immediately add to feature counts
        4. Delete the temporary sorted BAM before returning

    Peak memory is now proportional to the number of unique
    (feature_type, feature_name) pairs rather than the number of reads.
    """
    import pysam
    import os
    import subprocess
    import tempfile

    feature_counts = defaultdict(float)    # (ft, fn) -> fractional count
    trf_counts     = defaultdict(float)    # fn       -> tRNA_tRF count
    tirna_counts   = defaultdict(float)    # fn       -> tRNA_tiRNA count

    # ------------------------------------------------------------------
    # Step 1: Sort BAM by query name into a temp file
    # samtools sort -n is memory-efficient (uses disk for merge passes)
    # ------------------------------------------------------------------
    tmp_dir      = os.path.dirname(bam_path)
    sorted_bam   = bam_path.replace(".bam", "_namesorted.bam")

    print(f"  Sorting BAM by query name -> {sorted_bam}", file=sys.stderr)

    result = subprocess.run(
        ["samtools", "sort", "-n", "-@", "2",
         "-o", sorted_bam, bam_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"[ERROR] samtools sort failed: {result.stderr}", file=sys.stderr)
        raise RuntimeError("samtools sort -n failed")

    # ------------------------------------------------------------------
    # Step 2: Stream through name-sorted BAM one read group at a time
    # ------------------------------------------------------------------
    bam  = pysam.AlignmentFile(sorted_bam, "rb")
    refs = bam.references

    current_rname  = None
    current_hits   = set()       # (feature_type, feature_name) for current read
    current_length = 0
    current_abund  = 1

    def flush_read(rname, hits, length, abund):
        """Distribute a single read's abundance fractionally."""
        if not hits:
            return
        n_hits = len(hits)
        share  = abund / n_hits
        for (ft, fn) in hits:
            # miRNA length filter — only count reads 20-22 nt
            # Reads mapping to miRNA outside this window are silently skipped.
            # All other feature types are counted regardless of length.
            if ft == "miRNA" and not (20 <= length <= 22):
                continue
            feature_counts[(ft, fn)] += share
            if ft == "tRNA":
                if TRF_MIN <= length <= TRF_MAX:
                    trf_counts[fn] += share
                elif TIRNA_MIN <= length <= TIRNA_MAX:
                    tirna_counts[fn] += share

    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped:
            continue

        rname   = aln.query_name
        ref     = refs[aln.reference_id]
        hit_key = label_map.get(ref, ("unknown", ref))
        length  = aln.query_length
        abund   = get_collapsed_count(rname)

        if rname != current_rname:
            # Flush the previous read group
            if current_rname is not None:
                flush_read(current_rname, current_hits,
                           current_length, current_abund)
            # Start a new read group
            current_rname  = rname
            current_hits   = {hit_key}
            current_length = length
            current_abund  = abund
        else:
            current_hits.add(hit_key)

    # Flush the final read group
    if current_rname is not None:
        flush_read(current_rname, current_hits,
                   current_length, current_abund)

    bam.close()

    # ------------------------------------------------------------------
    # Step 3: Remove temporary sorted BAM
    # ------------------------------------------------------------------
    try:
        os.remove(sorted_bam)
        print(f"  Removed temporary sorted BAM", file=sys.stderr)
    except OSError as e:
        print(f"[WARN] Could not remove {sorted_bam}: {e}", file=sys.stderr)

    n_features = len(feature_counts)
    n_types    = len({ft for (ft, fn) in feature_counts})
    print(f"  Counted {n_features} feature-name entries across "
          f"{n_types} feature types", file=sys.stderr)

    return feature_counts, trf_counts, tirna_counts

# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def write_counts(sample_id, genome, tissue,
                 feature_counts, trf_counts, tirna_counts,
                 total_library_reads, out_path):
    """
    Write per-feature-name counts to TSV.

    Columns: library, genome, tissue, feature_type, feature_name,
             raw_count, RPM

    tRNA is expanded into three rows per sequence:
        tRNA_all, tRNA_tRF, tRNA_tiRNA
    All other feature types appear as-is.
    """
    rpm_scale = (1_000_000 / total_library_reads
                 if total_library_reads > 0 else 0.0)

    # Collect all tRNA feature names for the sub-category rows
    tRNA_names = sorted({fn for (ft, fn) in feature_counts if ft == "tRNA"})

    rows = []

    # Group by feature type, sort feature types alphabetically
    # Separate tRNA for special handling
    by_type = defaultdict(list)    # ft -> [(fn, count)]
    for (ft, fn), count in feature_counts.items():
        by_type[ft].append((fn, count))

    for ft in sorted(by_type.keys()):
        entries = sorted(by_type[ft], key=lambda x: x[0])

        if ft == "tRNA":
            # tRNA_all rows
            for fn, count in entries:
                rows.append(("tRNA_all", fn, count))

            # tRNA_tRF rows
            for fn in tRNA_names:
                count = trf_counts.get(fn, 0.0)
                rows.append(("tRNA_tRF", fn, count))

            # tRNA_tiRNA rows
            for fn in tRNA_names:
                count = tirna_counts.get(fn, 0.0)
                rows.append(("tRNA_tiRNA", fn, count))
        else:
            for fn, count in entries:
                rows.append((ft, fn, count))

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow([
            "library", "genome", "tissue",
            "feature_type", "feature_name",
            "raw_count", "RPM"
        ])
        for (ft, fn, count) in rows:
            writer.writerow([
                sample_id,
                genome,
                tissue,
                ft,
                fn,
                f"{count:.4f}",
                f"{count * rpm_scale:.4f}"
            ])

    print(f"[{sample_id}] Written {len(rows)} feature rows -> {out_path}",
          file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # CHANGED: pass genome_name from args rather than parsing from sample_id
    genome, tissue = parse_library_metadata(args.sample_id, args.genome_name)

    print(f"[{args.sample_id}] Genome: {genome}  Tissue: {tissue}",
          file=sys.stderr)

    # CHANGED: build label map from BAM header, not from FASTA file
    label_map           = build_label_map_from_bam(args.bam)
    label_map           = enrich_label_map_with_full_headers(label_map, args.combined_fa)
    total_library_reads = parse_bowtie2_total_reads(args.bowtie2_stats)

    print(f"[{args.sample_id}] Library reads: {total_library_reads}",
          file=sys.stderr)

    feature_counts, trf_counts, tirna_counts = \
        collect_feature_counts(args.bam, label_map)

    write_counts(
        args.sample_id, genome, tissue,
        feature_counts, trf_counts, tirna_counts,
        total_library_reads,
        args.out
    )


if __name__ == "__main__":
    main()