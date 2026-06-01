#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

// Modules
include { fastqc                        } from './modules/fastqc'
include { multiqc                       } from './modules/multiqc'
include { trimGalore                    } from './modules/trimGalore'
include { buildIndex                    } from './modules/buildIndex'
include { annotate_tRNA                 } from './modules/annotate_tRNA'
include { annotate_rRNA                 } from './modules/annotate_rRNA'
include { portAnnotations               } from './modules/portAnnotations'
include { align_sRNA                    } from './modules/align_sRNA'
include { extractSizeDistribution       } from './modules/extractSizeDistribution'
include { mergeSizeDistributions        } from './modules/mergeSizeDistributions'
include { featureCounts_sRNA            } from './modules/featureCounts'
include { quantify_sRNA_diversity       } from './modules/quantify_srna_diversity'
include { buildCombinedAnnotIndex       } from './modules/build_combined_annot_index'
include { alignToCombinedAnnotations    } from './modules/align_to_combined_annotations'
include { mergeDiversityOutputs         } from './modules/merge_diversity_outputs'
include { computeBetaDiversity          } from './modules/compute_beta_diversity'

// Parameters: see ./nextflow.config

// Pre-flight checks (Header)
log.info """\
 ===================================
 PIPELINE: Ex-sRNA-NF
 AUTHOR  : Scott Lewis
 ===================================
 reads            : ${params.reads}
 genomes          : ${params.genomes}
 outdir           : ${params.outdir}
 min_length       : ${params.min_length}
 max_length       : ${params.max_length}
 ===================================
 """

// Create channels
// ch_sRNA_reads = Channel.fromPath(params.reads, checkIfExists: true)
ch_genomes = Channel.fromPath(params.genomes, checkIfExists: true)
ch_sRNA_reads = Channel.fromPath(params.reads, checkIfExists: true)
    .filter { fastq ->
        def exclude = params.exclude_prefix ? params.exclude_prefix.split(",") : []
        !exclude.any { prefix -> fastq.name.startsWith(prefix) }
    }

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

def matchLibraryToIndex(library, indices_list) {
    try {
        // Get the base filename without path and extension
        def libraryBaseName = library.simpleName.toString()
        
        // Remove common suffixes added by trimming tools
        libraryBaseName = libraryBaseName.replaceAll(/_trimmed$/, '')
        libraryBaseName = libraryBaseName.replaceAll(/_val_[12]$/, '')
        
        // Split on underscore and get the first part as the genome prefix
        def libraryNameParts = libraryBaseName.split('_')
        
        if (libraryNameParts.size() == 0) {
            log.warn "Could not parse library name: ${library.name}"
            return null
        }
        
        def genome_prefix = libraryNameParts[0]
        
        // Log for debug: print the structure
        log.debug "Processing library: ${libraryBaseName}, looking for genome: ${genome_prefix}"
        log.debug "Indices list type: ${indices_list.getClass()}, size: ${indices_list.size()}"
        
        // indices_list is a list of tuples: [[genome_name1, [files1]], [genome_name2, [files2]], ...]
        def matching_index = null
        
        indices_list.eachWithIndex { index_item, idx ->
            log.debug "Index ${idx}: type=${index_item.getClass()}, content=${index_item}"
            
            // Handles both list and tuple formats
            def genome_name = index_item[0].toString()
            def index_files = index_item[1]
            
            log.debug "  Comparing '${genome_name}' with '${genome_prefix}'"
            
            if (genome_name == genome_prefix) {
                matching_index = [genome_name, index_files]
                log.debug "  MATCH FOUND!"
            }
        }
        
        if (matching_index == null) {
            def available = indices_list.collect { it[0].toString() }.join(', ')
            log.warn "No matching index found for library: ${libraryBaseName} (looking for genome prefix: ${genome_prefix})"
            log.warn "Available genome indices: ${available}"
            return null
        }
        
        log.info "Matched library ${libraryBaseName} to genome ${matching_index[0]}"
        return tuple(library, matching_index[0], matching_index[1])
        
    } catch (Exception e) {
        log.error "Error in matchLibraryToIndex: ${e.message}"
        log.error "Library: ${library}"
        log.error "Indices list: ${indices_list}"
        e.printStackTrace()
        return null
    }
}

// =============================================================================
// WORKFLOW
// =============================================================================
workflow {
    file(params.outdir).mkdirs()

    // --- QC ---
    ch_fastqc = fastqc(ch_sRNA_reads)
    multiqc(ch_fastqc.collect())

    // --- Trimming ---
    trimGalore(ch_sRNA_reads)
    ch_trimmed = trimGalore.out.trimmed_reads

    // --- Genome indexing ---
    ch_indices = buildIndex(ch_genomes)

    // --- Annotation ---
    annotate_tRNA(ch_genomes)
    annotate_rRNA(ch_genomes)
    portAnnotations(ch_genomes)

    // --- Combined annotation index ---
    ch_annot_ready = annotate_tRNA.out.tRNA_fasta
        .map { sample_id, fasta -> tuple(sample_id.replace('_genome', ''), fasta) }
        .mix(annotate_rRNA.out.rRNA_fasta
            .map { sample_id, fasta -> tuple(sample_id.replace('_genome', ''), fasta) })
//        .mix(portAnnotations.out.miRNA_fasta)
//        .mix(portAnnotations.out.TAS_fasta)
//        .mix(portAnnotations.out.TE_fasta)
//        .mix(portAnnotations.out.UTR5_fasta)
//        .mix(portAnnotations.out.UTR3_fasta)
//        .mix(portAnnotations.out.CDS_fasta)
//        .mix(portAnnotations.out.lncRNA_fasta)
//        .mix(portAnnotations.out.snRNA_fasta)
//        .mix(portAnnotations.out.snoRNA_fasta)
        .mix(
        // transpose splits tuple(sample_id, [fa1, fa2, fa3, ...])
        // into individual tuple(sample_id, fa1), tuple(sample_id, fa2), etc.
        portAnnotations.out.annotation_fastas
            .transpose()
        )
        .groupTuple(by: 0)

    buildCombinedAnnotIndex(ch_annot_ready)

    // --- Match trimmed reads to combined annotation index ---
    ch_trimmed
        .combine(buildCombinedAnnotIndex.out.index)
        .map { fastq, sample_id, combined_fa, index_files ->
            def libBase   = fastq.simpleName
                .replaceAll(/_trimmed$/, '')
                .replaceAll(/_val_[12]$/, '')
            def libPrefix = libBase.split('_')[0]
            if (sample_id.startsWith(libPrefix)) {
                log.info "Matched ${libBase} to combined annotation index ${sample_id}"
                return tuple(libBase, sample_id, fastq, combined_fa, index_files)
            }
            return null
        }
        .filter { it != null }
        .set { ch_matched_combined_annot }

    alignToCombinedAnnotations(ch_matched_combined_annot)

    // --- Genome alignment ---
    ch_trimmed
        .combine(ch_indices)
        .map { library, genome_name, index_files ->
            def libraryBaseName = library.simpleName.toString()
            libraryBaseName = libraryBaseName.replaceAll(/_trimmed$/, '')
            libraryBaseName = libraryBaseName.replaceAll(/_val_[12]$/, '')
            def libraryNameParts = libraryBaseName.split('_')
            def genome_prefix = libraryNameParts[0]
            if (genome_name.startsWith(genome_prefix)) {
                log.info "Matched library ${libraryBaseName} to genome ${genome_name}"
                return tuple(library, genome_name, index_files)
            } else {
                log.debug "Skipping: ${libraryBaseName} (prefix: ${genome_prefix}) doesn't match ${genome_name}"
                return null
            }
        }
        .filter { it != null }
        .set { ch_matched }

    align_sRNA(ch_matched)

    // --- Key genome BAMs for joining ---
    ch_genome_bams_keyed = align_sRNA.out.bam
        .map { bam ->
            def lib_name = bam.simpleName
                .replaceAll(/_vs_.*/, '')
                .replaceAll(/_trimmed$/, '')
            tuple(lib_name, bam)
        }

    ch_genome_bais_keyed = align_sRNA.out.index
        .map { bai ->
            def lib_name = bai.simpleName
                .replaceAll(/_vs_.*/, '')
                .replaceAll(/_trimmed$/, '')
                .replaceAll(/\.bam$/, '')
            tuple(lib_name, bai)
        }

    // --- Key annotation BAMs for joining ---
    ch_annot_bams_keyed = alignToCombinedAnnotations.out.bam
        .map { lib_name, bam, bai, combined_fa, stats ->
            tuple(lib_name, bam, bai)
        }

    // --- Join annotation + genome BAMs for quantification ---
    ch_for_quantify = alignToCombinedAnnotations.out.bam
        .map { lib_name, annot_bam, annot_bai, combined_fa, stats ->
            tuple(lib_name, annot_bam, annot_bai, combined_fa, stats)
        }
        .join(ch_genome_bams_keyed, by: 0)
        .join(ch_genome_bais_keyed, by: 0)
        .map { lib_name, annot_bam, annot_bai, combined_fa, stats,
               genome_bam, genome_bai ->
            tuple(lib_name, annot_bam, annot_bai, combined_fa, stats,
                  genome_bam, genome_bai)
        }

    // --- Quantification and diversity ---
    quantify_sRNA_diversity(ch_for_quantify)

    mergeDiversityOutputs(
        quantify_sRNA_diversity.out.diversity.collect(),
        quantify_sRNA_diversity.out.metrics.collect()
    )

    computeBetaDiversity(
        mergeDiversityOutputs.out.diversity
    )

    // --- Size distributions ---
    ch_for_size_dist = ch_genome_bams_keyed
        .join(ch_genome_bais_keyed, by: 0)
        .join(ch_annot_bams_keyed,  by: 0)
        .map { lib_name, genome_bam, genome_bai, annot_bam, annot_bai ->
            tuple(genome_bam, genome_bai, annot_bam, annot_bai)
        }

    extractSizeDistribution(ch_for_size_dist)

    mergeSizeDistributions(
        extractSizeDistribution.out.size_dist_counts.collect(),
        extractSizeDistribution.out.size_dist_rpm.collect()
    )
}

workflow.onComplete {
    log.info """\
    ===================================
    Pipeline completed!
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Results: ${params.outdir}
    ===================================
    """
}