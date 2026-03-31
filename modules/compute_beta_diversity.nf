process computeBetaDiversity {
    cpus 1
    memory '4 GB'
    publishDir { "${params.outdir}/06_srna_diversity" }, mode: 'copy', overwrite: false

    input:
    path diversity_tsv    // i.e. 'all_samples_srna_diversity.tsv'

    output:
    path "bray_curtis_matrix.tsv",    emit: bray_curtis
    path "aitchison_matrix.tsv",      emit: aitchison
    path "beta_diversity_summary.tsv", emit: summary

    script:
    """
    compute_beta_diversity.py \
        --diversity     ${diversity_tsv} \
        --out_bray_curtis  bray_curtis_matrix.tsv \
        --out_aitchison    aitchison_matrix.tsv \
        --out_summary      beta_diversity_summary.tsv
    """
}