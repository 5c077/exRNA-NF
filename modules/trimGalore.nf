process trimGalore {
    tag "${reads.simpleName}"
    cpus 4
    memory '4 GB'
    publishDir "${params.outdir}/02_trimGalore", mode: 'symlink', overwrite: true

    input:
    path reads

    output:
    path "*_trimmed.fq.gz", emit: trimmed_reads
    path "*_trimming_report.txt", emit: reports
    path "*_fastqc.{zip,html}", emit: fastqc

    script:
    def adapter_cmd = ""
    if (params.use_small_rna_adapter) {
        adapter_cmd = "--small_rna"
    } else if (params.adapter) {
        adapter_cmd = "${params.adapter}"
    }
    
    """
    trim_galore \
        --fastqc \
        --fastqc_args "--outdir ./" \
        -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC \
        --clip_R1 1 \
        --length ${params.min_length} \
        --max_length ${params.max_length} \
        --max_n 0 \
        --cores ${task.cpus} \
        ${reads}
    """
}