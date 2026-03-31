process extractSizeDistribution {
    tag "${genome_bam.simpleName}"
    cpus 2
    memory '4 GB'
    publishDir "${params.outdir}/04_size_distribution", mode: 'symlink', overwrite: true

    input:
    // genome_bam    : BAM from align_sRNA (all reads mapped to genome)
    // annot_bam     : BAM from ALIGN_TO_COMBINED_ANNOTATIONS (reads mapped to
    //                 combined feature index, reference names = featureType|seqId)
    tuple path(genome_bam), path(genome_bai),
          path(annot_bam),  path(annot_bai)

    output:
    path "${genome_bam.simpleName}_size_dist_counts.csv", emit: size_dist_counts
    path "${genome_bam.simpleName}_size_dist_rpm.csv",    emit: size_dist_rpm

    script:
    """
    #!/usr/bin/env python3
    import pysam
    import csv
    from collections import defaultdict

    genome_bam_path = "${genome_bam}"
    annot_bam_path  = "${annot_bam}"
    library_name    = "${genome_bam.simpleName}"

    # ------------------------------------------------------------------
    # Step 1: Total size distribution from genome BAM
    # ------------------------------------------------------------------
    total_counts = defaultdict(int)

    try:
        bam = pysam.AlignmentFile(genome_bam_path, "rb", check_sq=False)
        for read in bam:
            if not read.is_unmapped:
                total_counts[read.query_length] += 1
        bam.close()
    except Exception as e:
        print(f"Error reading genome BAM: {e}")

    total_reads = sum(total_counts.values())

    # ------------------------------------------------------------------
    # Step 2: Per-feature size distributions from combined annotation BAM
    # Reference names follow the format: featureType|sequenceId
    # A read may align to multiple feature types (fractional counting);
    # for size distribution purposes we count it once per feature type hit.
    # ------------------------------------------------------------------
    feature_counts = defaultdict(lambda: defaultdict(int))   # {feature: {length: count}}
    read_features  = defaultdict(set)                         # {read_name: {feature_types}}
    read_lengths   = {}                                       # {read_name: length}

    try:
        bam = pysam.AlignmentFile(annot_bam_path, "rb", check_sq=False)
        refs = bam.references
        for read in bam:
            if read.is_unmapped:
                continue
            ref          = refs[read.reference_id]
            feature_type = ref.split("|")[0]
            rname        = read.query_name
            read_features[rname].add(feature_type)
            read_lengths[rname] = read.query_length
        bam.close()
    except Exception as e:
        print(f"Error reading annotation BAM: {e}")

    # Assign each read to every feature type it mapped to
    # (mirrors the fractional counting logic in quantify_srna_diversity.py)
    for rname, features in read_features.items():
        length = read_lengths[rname]
        for ft in features:
            feature_counts[ft][length] += 1

    # ------------------------------------------------------------------
    # Step 3: Determine global size range across all categories
    # ------------------------------------------------------------------
    all_lengths = set(total_counts.keys())
    for ft_counts in feature_counts.values():
        all_lengths.update(ft_counts.keys())

    if not all_lengths:
        # Write empty files so the pipeline does not fail
        for suffix in ["_size_dist_counts.csv", "_size_dist_rpm.csv"]:
            with open(f"{library_name}{suffix}", "w") as f:
                f.write("library,Category\\n")
        raise SystemExit(0)

    min_size   = min(all_lengths)
    max_size   = max(all_lengths)
    size_range = list(range(min_size, max_size + 1))

    # ------------------------------------------------------------------
    # Step 4: Build all rows — Total first, then feature types alphabetically
    # ------------------------------------------------------------------
    # Each row: (library_name, category_label, {length: count})
    rows = [("Total", total_counts)]
    for ft in sorted(feature_counts.keys()):
        rows.append((ft, feature_counts[ft]))

    header = ["library", "Category"] + [str(s) for s in size_range]

    # Write counts CSV
    with open(f"{library_name}_size_dist_counts.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for category, counts in rows:
            row = [library_name, category]
            row += [counts.get(s, 0) for s in size_range]
            writer.writerow(row)

    # Write RPM CSV (denominator = total reads from genome BAM)
    with open(f"{library_name}_size_dist_rpm.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for category, counts in rows:
            row = [library_name, category]
            if total_reads > 0:
                row += [round((counts.get(s, 0) / total_reads) * 1_000_000, 2)
                        for s in size_range]
            else:
                row += [0] * len(size_range)
            writer.writerow(row)
    """
}