process portAnnotations {
    tag "${genome_fasta.simpleName}"
    publishDir "${params.outdir}/annotations/hairpin", mode: 'copy', overwrite: true, pattern: "*_hairpin.fa"
    publishDir "${params.outdir}/annotations/miRNA",   mode: 'copy', overwrite: true, pattern: "*_miRNA.fa"
    publishDir "${params.outdir}/annotations/TAS",     mode: 'copy', overwrite: true, pattern: "*_TAS.fa"
    publishDir "${params.outdir}/annotations/TE",      mode: 'copy', overwrite: true, pattern: "*_TE.fa"
    publishDir "${params.outdir}/annotations/cDNA",    mode: 'copy', overwrite: true, pattern: "*_cDNA.fa"

    input:
    path genome_fasta

    output:
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_hairpin.fa"), emit: hairpin_fasta, optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_miRNA.fa"),   emit: miRNA_fasta,   optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_TAS.fa"),     emit: TAS_fasta,     optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_TE.fa"),      emit: TE_fasta,      optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_cDNA.fa"),    emit: cDNA_fasta,    optional: true

    script:
    def base = genome_fasta.simpleName.replace('_genome', '')
    def genome_dir = genome_fasta.toRealPath().parent
    """
    echo "DEBUG: genome_fasta = ${genome_fasta}"
    echo "DEBUG: base         = ${base}"
    echo "DEBUG: genome_dir   = ${genome_dir}"
    echo "DEBUG: resolved parent contents:"
    ls -la "${genome_dir}/" || echo "Cannot list genome_dir"

    for suffix in _hairpin.fa _miRNA.fa _TAS.fa _TE.fa _cDNA.fa; do
        src="${genome_dir}/${base}\${suffix}"
        echo "DEBUG: looking for \${src}"
        if [ -f "\${src}" ]; then
            echo "DEBUG: FOUND \${src}, copying..."
            cp "\${src}" "${base}\${suffix}"
        else
            echo "DEBUG: NOT FOUND \${src}"
        fi
    done
    """
}