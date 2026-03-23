process annotate_tRNA {
    tag "${genome_fasta.simpleName}"
    cpus 4
    memory '8 GB'
    publishDir "${params.outdir}/annotations/tRNA", mode: 'symlink', overwrite: false

    input:
    path genome_fasta

    output:
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}_tRNA.out"), emit: tRNA_table
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}_tRNA.bed"), emit: tRNA_bed
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}_tRNA.fa"), emit: tRNA_fasta
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}_tRNA.ss"), emit: tRNA_structure

    script:
    """
    # Run tRNAscan-SE
    # -o: output table
    # -f: output secondary structure
    # -b: output BED format
    # -a: output fasta
    # --thread: number of threads
    tRNAscan-SE \
        -o ${genome_fasta.simpleName}_tRNA.out \
        -f ${genome_fasta.simpleName}_tRNA.ss \
        -b ${genome_fasta.simpleName}_tRNA.bed \
        -a ${genome_fasta.simpleName}_tRNA.fa \
        --thread ${task.cpus} \
        ${genome_fasta}
    """
}