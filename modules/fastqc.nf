process fastqc {
    tag "${reads.simpleName}"
    cpus 4
    memory '4 GB'
    publishDir "${params.outdir}/00_fastqc", mode: 'symlink', overwrite: false

    input:
    path reads

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}