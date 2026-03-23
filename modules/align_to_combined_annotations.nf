process alignToCombinedAnnotations {
    tag "${lib_name}"
    cpus 8
    memory '24 GB'
    publishDir { "${params.outdir}/05_annotation_alignment/combined" }, mode: 'symlink', overwrite: false

    input:
    tuple val(lib_name), val(sample_id), path(fastq), path(combined_fa), path(index_dir)

    output:
    tuple val(lib_name),
      path("${lib_name}_vs_combined_annotations.bam"),
      path("${lib_name}_vs_combined_annotations.bam.csi"),
      path(combined_fa),                             
      path("${lib_name}_bowtie2_stats.txt"),         
      emit: bam

    script:
    def index_prefix = "${index_dir}/${sample_id}_combined"
    """
    bowtie2 \
        -q \
        -N 0 \
        -k ${params.max_feature_types ?: 10} \
        --no-1mm-upfront \
        --no-unal \
        -p ${task.cpus} \
        -x ${index_prefix} \
        -U ${fastq} \
        2> ${lib_name}_bowtie2_stats.txt \
    | samtools sort \
        -@ ${task.cpus} \
        -o ${lib_name}_vs_combined_annotations.bam

    if [ ! -s "${lib_name}_vs_combined_annotations.bam" ]; then
        echo "ERROR: bowtie2 produced no output for ${lib_name}"
        echo "bowtie2 log:"
        cat ${lib_name}_bowtie2_stats.txt
        exit 1
    fi

    samtools index \
        -@ ${task.cpus} \
        -c ${lib_name}_vs_combined_annotations.bam

    echo "bowtie2 alignment summary for ${lib_name}:"
    cat ${lib_name}_bowtie2_stats.txt
    """
}