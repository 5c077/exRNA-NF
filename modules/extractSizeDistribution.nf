process extractSizeDistribution {
    tag "${genome_bam.simpleName}"
    cpus 2
    memory '8 GB'
    maxForks 16            // ADDED: cap at 16 concurrent jobs (16 * 8 GB = 128 GB peak)
    publishDir { "${params.outdir}/04_size_distribution" }, mode: 'symlink', overwrite: false

    input:
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
    import re
    from collections import defaultdict

    genome_bam_path = "${genome_bam}"
    annot_bam_path  = "${annot_bam}"
    library_name    = "${genome_bam.simpleName}"

    # ------------------------------------------------------------------
    # Parse genome and tissue type from library name
    # Library name pattern:
    #   {prefix}_trimmed_vs_{genome}_genome_final
    # e.g. Col0_Col0_CL_R2_S28_R1_001_trimmed_vs_Col0_TAIR10_genome_final
    # ------------------------------------------------------------------

    # Extract genome: everything between _vs_ and _genome_final
    genome_match = re.search(r'_vs_(.+?)_genome_final', library_name)
    genome       = genome_match.group(1) if genome_match else "unknown"

    # Extract tissue type by scanning for known treatment tokens
    # in the portion of the name before _vs_
    lib_prefix   = library_name.split("_vs_")[0]
    tissue_tokens = ["AWF", "CL", "LSW", "LSF"]
    tissue = "unknown"
    for token in tissue_tokens:
        # Match token as a whole underscore-delimited field
        if re.search(r'(?:^|_)' + token + r'(?:_|\$)', lib_prefix):
            tissue = token
            break

    print(f"Library: {library_name}")
    print(f"  Genome:  {genome}")
    print(f"  Tissue:  {tissue}")

    # ------------------------------------------------------------------
    # Step 1: Total size distribution from genome BAM
    # CHANGED: also collect genome read names in same pass for Other computation
    # ------------------------------------------------------------------
    total_counts      = defaultdict(int)
    genome_read_names = set()                    # ADDED: collect names here

    try:
        bam = pysam.AlignmentFile(genome_bam_path, "rb", check_sq=False)
        for read in bam:
            if not read.is_unmapped:
                total_counts[read.query_length] += 1
                genome_read_names.add(read.query_name)   # ADDED
        bam.close()
    except Exception as e:
        print(f"Error reading genome BAM: {e}")

    total_reads = sum(total_counts.values())

    # ------------------------------------------------------------------
    # Step 2: Per-feature size distributions from combined annotation BAM
    # UNCHANGED
    # ------------------------------------------------------------------
    feature_counts = defaultdict(lambda: defaultdict(int))
    read_features  = defaultdict(set)
    read_lengths   = {}

    try:
        bam  = pysam.AlignmentFile(annot_bam_path, "rb", check_sq=False)
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

    for rname, features in read_features.items():
        length = read_lengths[rname]
        for ft in features:
            feature_counts[ft][length] += 1

    # ------------------------------------------------------------------
    # Step 3: Other — derive from set difference, no second BAM pass
    # CHANGED: use genome_read_names set instead of re-reading the BAM
    # ------------------------------------------------------------------
    other_read_names = genome_read_names - set(read_features.keys())  # CHANGED

    if other_read_names:
        # We need lengths for Other reads — collect from read_lengths where available,
        # but Other reads by definition are NOT in read_lengths (annotation BAM only).
        # So we need a targeted second pass only for Other read lengths,
        # or use the total_counts distribution as a proxy.
        # 
        # Most accurate: single targeted pass collecting lengths for unmapped-to-annot reads
        other_counts = defaultdict(int)
        try:
            bam = pysam.AlignmentFile(genome_bam_path, "rb", check_sq=False)
            for read in bam:
                if not read.is_unmapped and read.query_name in other_read_names:
                    other_counts[read.query_length] += 1
            bam.close()
        except Exception as e:
            print(f"Error computing Other lengths: {e}")

        if other_counts:
            feature_counts["Other"] = other_counts

    # Free memory before writing output
    del genome_read_names
    del read_features
    del read_lengths

    # ------------------------------------------------------------------
    # Step 4: Global size range
    # ------------------------------------------------------------------
    all_lengths = set(total_counts.keys())
    for ft_counts in feature_counts.values():
        all_lengths.update(ft_counts.keys())

    if not all_lengths:
        for suffix in ["_size_dist_counts.csv", "_size_dist_rpm.csv"]:
            with open(f"{library_name}{suffix}", "w") as f:
                f.write("library,genome,tissue,feature_type\\n")
        raise SystemExit(0)

    min_size   = min(all_lengths)
    max_size   = max(all_lengths)
    size_range = list(range(min_size, max_size + 1))

    # ------------------------------------------------------------------
    # Step 5: Build rows
    # Order: Total, Other, then feature types alphabetically
    # CHANGED: rows now include genome and tissue columns
    # ------------------------------------------------------------------
    rows = [("Total", total_counts)]

    if "Other" in feature_counts:
        rows.append(("Other", feature_counts.pop("Other")))

    for ft in sorted(feature_counts.keys()):
        rows.append((ft, feature_counts[ft]))

    # CHANGED: header now includes genome and tissue columns
    header = ["library", "genome", "tissue", "feature_type"] + \
             [str(s) for s in size_range]

    # ------------------------------------------------------------------
    # Step 6: Write counts CSV
    # CHANGED: rows now include genome and tissue values
    # ------------------------------------------------------------------
    with open(f"{library_name}_size_dist_counts.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for feature_type, counts in rows:
            row = [library_name, genome, tissue, feature_type]
            row += [counts.get(s, 0) for s in size_range]
            writer.writerow(row)

    # ------------------------------------------------------------------
    # Step 7: Write RPM CSV
    # CHANGED: rows now include genome and tissue values
    # ------------------------------------------------------------------
    with open(f"{library_name}_size_dist_rpm.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for feature_type, counts in rows:
            row = [library_name, genome, tissue, feature_type]
            if total_reads > 0:
                row += [round((counts.get(s, 0) / total_reads) * 1_000_000, 2)
                        for s in size_range]
            else:
                row += [0] * len(size_range)
            writer.writerow(row)
    """
}