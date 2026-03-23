process buildIndex {
    tag "${genome_fasta.simpleName}"
    cpus 8
    memory '8 GB'
    publishDir "${params.outdir}/genome_indices", mode: 'symlink', overwrite: false

    input:
    path genome_fasta

    output:
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}*.bt2*"), emit: indices

    script:
    """
    bowtie2-build --threads ${task.cpus} ${genome_fasta} ${genome_fasta.simpleName}
    """
}