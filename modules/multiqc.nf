process multiqc {
    cpus 2
    memory '2 GB'
    publishDir "${params.outdir}/01_multiqc", mode: 'symlink', overwrite: false
    
    input:
    path 'fastqc/*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc fastqc/
    """
}