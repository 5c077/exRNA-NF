process annotate_rRNA {
    tag "${genome_fasta.simpleName}"
    cpus 4
    memory '8 GB'
    publishDir "${params.outdir}/annotations/rRNA", mode: 'symlink', overwrite: true
    
    input:
    path genome_fasta
    
    output:
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}_rRNA.gff"), emit: rRNA_gff
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}_rRNA.fa"), emit: rRNA_fasta
    tuple val("${genome_fasta.simpleName}"), path("${genome_fasta.simpleName}_rRNA.xml"), emit: rRNA_xml, optional: true
    
    script:
    """
    # Detect kingdom based on genome name (default to eukaryote)
    KINGDOM="euk"
    if [[ "${genome_fasta.simpleName}" =~ [Bb]ac ]]; then
        KINGDOM="bac"
    elif [[ "${genome_fasta.simpleName}" =~ [Aa]rc ]]; then
        KINGDOM="arc"
    fi
    
    # Simplify FASTA headers (remove everything after first space)
    # This prevents issues with nhmmer alphabet detection
    sed 's/ .*//' ${genome_fasta} > ${genome_fasta.simpleName}_clean.fa
    
    # Run barrnap with --incseq to embed sequences in GFF3
    barrnap \
        --kingdom \${KINGDOM} \
        --threads ${task.cpus} \
        --incseq \
        --quiet \
        --outseq ${genome_fasta.simpleName}_rRNA.fa \
        ${genome_fasta.simpleName}_clean.fa \
        > ${genome_fasta.simpleName}_rRNA_with_seq.gff || true
    
    echo "Using Kingdom: \${KINGDOM}"
    
    # Extract GFF annotations (everything before ##FASTA)
    awk '/^##FASTA/,0{exit} {print}' ${genome_fasta.simpleName}_rRNA_with_seq.gff > ${genome_fasta.simpleName}_rRNA.gff 2>/dev/null || true
    
    # Ensure GFF and FASTA files exist even if no rRNAs found
    touch ${genome_fasta.simpleName}_rRNA.gff
    touch ${genome_fasta.simpleName}_rRNA.fa
    
    # Clean up intermediate files
    #rm -f ${genome_fasta.simpleName}_clean.fa
    #rm -f ${genome_fasta.simpleName}_rRNA_with_seq.gff
    """
}