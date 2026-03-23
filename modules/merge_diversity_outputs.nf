process mergeDiversityOutputs {
    cpus 1
    memory '4 GB'
    publishDir { "${params.outdir}/06_srna_diversity" }, mode: 'copy', overwrite: true

    input:
    path diversity_files
    path metrics_files

    output:
    path "all_samples_srna_diversity.tsv",    emit: diversity
    path "all_samples_alpha_diversity.tsv",   emit: alpha
    path "all_samples_beta_diversity.tsv",    emit: beta

    script:
    def w_seg_arg = params.w_seg ? "--w_seg \"${params.w_seg}\"" : ""
    """
    merge_diversity_outputs.py \
        --diversity ${diversity_files} \
        --metrics   ${metrics_files} \
        --out_diversity all_samples_srna_diversity.tsv \
        --out_alpha     all_samples_alpha_diversity.tsv \
        --out_beta      all_samples_beta_diversity.tsv \
        ${w_seg_arg}
    """
}