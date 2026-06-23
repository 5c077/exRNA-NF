process COUNT_SRNA_FEATURES {
    tag "${lib_name}"
    cpus 4          // CHANGED: was 2; samtools sort -n uses -@ 2 threads
    memory '24 GB'  // CHANGED: was 8 GB; sort needs disk+RAM headroom
    maxForks 6      // ADDED: 6 * 24 GB = 144 GB peak, within 256 GB
    publishDir { "${params.outdir}/07_feature_counts" }, mode: 'copy', overwrite: true

    input:
    tuple val(lib_name),
          path(annot_bam),
          path(annot_bai),
          path(combined_fa),
          path(bowtie2_stats)

    output:
    path "${lib_name}_feature_counts.csv", emit: feature_counts

    script:
    // CHANGED: derive genome name from combined_fa filename
    // e.g. Lsa_Salinas_V15_Gapless_Final_combined.fa -> Lsa_Salinas_V15_Gapless_Final
    def genome_name = combined_fa.simpleName.replaceAll(/_combined$/, '')
    """
    count_srna_features.py \
        --sample_id     ${lib_name} \
        --genome_name   ${genome_name} \
        --bam           ${annot_bam} \
        --combined_fa   ${combined_fa} \
        --bowtie2_stats ${bowtie2_stats} \
        --out           ${lib_name}_feature_counts.csv
    """
}