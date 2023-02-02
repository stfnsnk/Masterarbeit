#!/usr/bin/env nextflow

nextflow.enable.dsl=2


//--------------------------------------------------------------------------------------//
//--------------------------------INCLUDING MODULES-------------------------------------//
//--------------------------------------------------------------------------------------//

include { fastp_paired; fastp_singleend }                           from "./modules/preprocessing.nf"
include { sam_to_sorted_bam; sam_to_sorted_filtered_bam }           from "./modules/samtools.nf"
include { hisat2_paired; minimap2; hisat2_singleend }               from "./modules/alignment.nf"
include { featureCounts_paired; featureCounts_long; featureCounts}  from "./modules/featureCounts.nf" 
include { multiqc; pycoQC }                                         from "./modules/quality_control.nf"
include { deseq2 }                                                  from "./modules/DESEQ2/DESEQ2.nf"






//--------------------------------------------------------------------------------------//
//------------------------------MAIN WORKFLOW SECTION-----------------------------------//
//--------------------------------------------------------------------------------------//


workflow ont_pipeline {
  log.info """\
        =============================================
        ============DGE Analysis Pipeline============
        ============== ONT LONG-READ ================
        =============================================
        Author: Stefan Senk - BIOINF22
        reads: ${params.reads_folder}
        reference genome: ${params.minimap2_index_folder}
        features (GTF): ${params.features_to_count_folder}
        sample info file: ${params.sample_info}
        """
        .stripIndent()

  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
    .view()
    .set{reads_ch}

  Channel
    .fromPath( "${params.minimap2_index_folder}/*.fa*" )
    .collect()
    .set { minimap_ref_ch } 

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }


  minimap2(reads_ch, minimap_ref_ch)

  sam_to_sorted_filtered_bam(minimap2.out.minimap2_sam, minimap_ref_ch)

  //pycoQC("${params.pyco_seq_summary}", "${params.pyco_barcode_summary}", sam_to_sorted_filtered_bam.out.out_bamfiles.collect())

  featureCounts_long(sam_to_sorted_filtered_bam.out.out_bamfiles, gene_annotation_file_ch)

  multiqc(featureCounts_long.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_long.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}

workflow short_read_pipeline_paired {
  log.info """\
         =============================================
         ============DGE Analysis Pipeline============
         ============PAIRED END SHORT-READ============
         =============================================
         Author: Stefan Senk - BIOINF22
         reads: ${params.reads_folder}
         reference genome: ${params.hisat2_index_folder}
         features (GTF): ${params.features_to_count_folder}
         sample info file: ${params.sample_info}
         """
         .stripIndent()

  Channel
    .fromFilePairs( "${params.reads_folder}/*_{R1,R2}*f*q.gz")
    .ifEmpty { exit 1, "no fastq files found at given path" }
    .set{ read_pair_ch }
  
  Channel
    .fromPath( "${params.hisat2_index_folder}/*.ht2" )
    .collect() 
    .map{ tuple( it.first().simpleName, it)}
    .ifEmpty { exit 1, "no hisat2 index files found - path was empty" }
    .set { hisat2_index_ch }  

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }


  fastp_paired(read_pair_ch)

  hisat2_paired(fastp_paired.out.trimmed_fastq_pair, hisat2_index_ch)
  
  sam_to_sorted_bam(hisat2_paired.out.sam_out_ch)

  featureCounts_paired(sam_to_sorted_bam.out.sorted_bamfiles, gene_annotation_file_ch)
 
  multiqc(featureCounts_paired.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_paired.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}

workflow short_read_pipeline_single {
  log.info """\
         =============================================
         ============DGE Analysis Pipeline============
         ============SINGLE END SHORT-READ============
         =============================================
         Author: Stefan Senk - BIOINF22
         reads: ${params.reads_folder}
         reference genome: ${params.hisat2_index_folder}
         features (GTF): ${params.features_to_count_folder}
         sample info file: ${params.sample_info}
         """
         .stripIndent()

  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file)}
    .view()
    .set{reads_ch}
  
  Channel
    .fromPath( "${params.hisat2_index_folder}/*.ht2" )
    .collect() 
    .map{ tuple( it.first().simpleName, it)}
    .ifEmpty { exit 1, "no hisat2 index files found - path was empty" }
    .set { hisat2_index_ch }

  fastp_singleend(reads_ch)

  hisat2_singleend(fastp_singleend.out.trimmed_fastq, hisat2_index_ch)
  
  sam_to_sorted_bam(hisat2_singleend.out.sam_out_ch)

  featureCounts(sam_to_sorted_bam.out.sorted_bamfiles, "${params.features_to_count_folder}/*.gtf*")

  multiqc(featureCounts.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}



