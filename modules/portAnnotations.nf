process portAnnotations {
    tag "${genome_fasta.simpleName}"
    publishDir "${params.outdir}/annotations/miRNA",   mode: 'copy', overwrite: true, pattern: "*_miRNA.fa"
    publishDir "${params.outdir}/annotations/TAS",     mode: 'copy', overwrite: true, pattern: "*_TAS.fa"
    publishDir "${params.outdir}/annotations/TE",      mode: 'copy', overwrite: true, pattern: "*_TE.fa"
    publishDir "${params.outdir}/annotations/5UTR",    mode: 'copy', overwrite: true, pattern: "*_5UTR.fa"
    publishDir "${params.outdir}/annotations/3UTR",    mode: 'copy', overwrite: true, pattern: "*_3UTR.fa"
    publishDir "${params.outdir}/annotations/CDS",     mode: 'copy', overwrite: true, pattern: "*_CDS.fa"
    publishDir "${params.outdir}/annotations/lncRNA",  mode: 'copy', overwrite: true, pattern: "*_lncRNA.fa"
    publishDir "${params.outdir}/annotations/snRNA",   mode: 'copy', overwrite: true, pattern: "*_snRNA.fa"
    publishDir "${params.outdir}/annotations/snoRNA",  mode: 'copy', overwrite: true, pattern: "*_snoRNA.fa"

    input:
    path genome_fasta

    output:
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_miRNA.fa"),   emit: miRNA_fasta,   optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_TAS.fa"),     emit: TAS_fasta,     optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_TE.fa"),      emit: TE_fasta,      optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_5UTR.fa"),    emit: UTR5_fasta,    optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_3UTR.fa"),    emit: UTR3_fasta,    optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_CDS.fa"),     emit: CDS_fasta,     optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_lncRNA.fa"),  emit: lncRNA_fasta,  optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_snRNA.fa"),   emit: snRNA_fasta,   optional: true
    tuple val("${genome_fasta.simpleName.replace('_genome', '')}"), path("*_snoRNA.fa"),  emit: snoRNA_fasta,  optional: true

    script:
    def base = genome_fasta.simpleName.replace('_genome', '')
    def genome_dir = genome_fasta.toRealPath().parent
    """
    echo "DEBUG: genome_fasta = ${genome_fasta}"
    echo "DEBUG: base         = ${base}"
    echo "DEBUG: genome_dir   = ${genome_dir}"
    echo "DEBUG: resolved parent contents:"
    ls -la "${genome_dir}/" || echo "Cannot list genome_dir"

    for suffix in _miRNA.fa _TAS.fa _TE.fa _5UTR.fa _3UTR.fa _CDS.fa _lncRNA.fa _snRNA.fa _snoRNA.fa; do
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