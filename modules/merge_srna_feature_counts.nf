process MERGE_SRNA_FEATURE_COUNTS {
    cpus 2
    memory '8 GB'
    publishDir { "${params.outdir}/07_feature_counts" }, mode: 'copy', overwrite: true

    input:
    path feature_count_files

    output:
    path "all_samples_feature_counts.csv", emit: merged

    script:
    """
    # Write header from the first file
    head -1 \$(ls *_feature_counts.csv | head -1) \
        > all_samples_feature_counts.csv

    # Concatenate all files, skipping the header line from each
    # awk NR>1 skips line 1 (the header) of each file
    # sort -k1,1 -k4,4 -k5,5 groups by library, feature_type, feature_name
    # for easier downstream analysis
    for f in *_feature_counts.csv; do
        tail -n +2 "\$f"
    done | sort -k1,1 -k4,4 -k5,5 -t ',' >> all_samples_feature_counts.csv

    # Report row count
    total=\$(wc -l < all_samples_feature_counts.csv)
    echo "Merged output: \${total} rows (including header)" >&2
    """
}