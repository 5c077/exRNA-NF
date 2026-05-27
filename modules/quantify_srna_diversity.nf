process quantify_sRNA_diversity {
    tag "${lib_name}"
    cpus 2
    memory '8 GB'
    publishDir { "${params.outdir}/06_srna_diversity" }, mode: 'copy', overwrite: true

    input:
    tuple val(lib_name),
          path(annot_bam),
          path(annot_bai),
          path(combined_fa),
          path(bowtie2_stats),
          path(genome_bam),       //Added for "other" calc.
          path(genome_bai)        //


    output:
    path "${lib_name}_srna_diversity.tsv",    emit: diversity
    path "${lib_name}_diversity_metrics.tsv", emit: metrics

    script:
    """
    quantify_srna_diversity.py \
        --sample_id     ${lib_name} \
        --bam           ${annot_bam} \
        --combined_fa   ${combined_fa} \
        --bowtie2_stats ${bowtie2_stats} \
        --genome_bam    ${genome_bam} \
        --out           ${lib_name}_srna_diversity.tsv \
        --metrics       ${lib_name}_diversity_metrics.tsv
    """
}