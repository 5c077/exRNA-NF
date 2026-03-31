process align_sRNA {
    tag "${sample_name}"
    cpus 8
    memory '24 GB'
    publishDir "${params.outdir}/03_alignment", mode: 'symlink', overwrite: true

    input:
    tuple path(fastq), val(genome_name), path(index_files)

    output:
    path "${sample_name}_final.bam", emit: bam
    path "${sample_name}_final.bam.csi", emit: index
    path "${sample_name}_alignment_stats.txt", emit: stats

    script:
    sample_name = "${fastq.simpleName}_vs_${genome_name}"
    """
    # Alignment with bowtie2
    bowtie2 \
        -q \
        -N 0 \
        -k 1 \
        --no-1mm-upfront \
        --no-unal \
        -x ${genome_name} \
        -U ${fastq} \
        -p ${task.cpus} \
        -S ${sample_name}.sam \
        2> ${sample_name}_alignment_stats.txt

    # Convert to BAM and sort
    #^Use this and uncomment below commands for more filtering
    samtools view -@ ${task.cpus} -bS ${sample_name}.sam | \
        samtools sort -@ ${task.cpus} -o ${sample_name}_final.bam

    # Mark and remove duplicates
    #samtools markdup -@ ${task.cpus} -r ${sample_name}_sorted.bam ${sample_name}_dedup.bam

    # Filter for properly mapped reads (forward and reverse strand)
    #samtools view -@ ${task.cpus} -h ${sample_name}_dedup.bam | \
    #    awk '\$0 ~ /^@/ || \$2 == 0 || \$2 == 16' | \
    #    samtools view -@ ${task.cpus} -bS - > ${sample_name}_final.bam

    # Index the final BAM using CSI format (handles large genomes)
    samtools index -@ ${task.cpus} -c ${sample_name}_final.bam

    # Clean up intermediate files
    #rm ${sample_name}.sam ${sample_name}_sorted.bam ${sample_name}_dedup.bam
    """
}