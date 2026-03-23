process featureCounts_sRNA {
    tag "${bam.simpleName}"
    cpus 4
    memory '8 GB'
    publishDir "${params.outdir}/05_counts", mode: 'symlink', overwrite: false

    input:
    path bam
    path annotation

    output:
    path "${bam.simpleName}_counts.txt", emit: counts
    path "${bam.simpleName}_counts.txt.summary", emit: summary

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -t miRNA \
        -g Name \
        -a ${annotation} \
        -o ${bam.simpleName}_counts.txt \
        ${bam}
    """
}