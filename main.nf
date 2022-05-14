#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.reads = "$baseDir/*.fastq.gz"
//params.hisat2_indexfiles ="$baseDir/ref_genome_grcm38/grcm38_snptran"

log.info """\
         reads: ${params.reads}
         hisat2_index: ${params.hisat2_index}
         """
         .stripIndent()


process HISAT2 {

  tag "$reads.baseName"
  
  container "${params.container_hisat}"

  input:
  path(reads)

  output:
  path('*.sam')

  script:
  """
  hisat2 -x "${params.hisat2_index}" \
         -U $reads \
         > ${reads.baseName}.sam
  """

}

workflow {

  Channel
    .fromPath( '/home/stefan/Dokumente/BIOINF22/Masterarbeit/Github/HISAT_versuch/*.fastq.gz' )
    .ifEmpty { exit 1, "reads was empty - no input files supplied" }
    .set { reads_ch }

  //CHannel für Index Files überhaupt nötig?
  // Channel.fromPath( '/home/stefan/Dokumente/BIOINF22/Masterarbeit/Github/HISAT_versuch/genome_snp_tran/*.ht2' )
  //   .ifEmpty{ exit 1, "no hisat2 index was found - no input files supplied" }
  //   .set( hisat2_index_ch )

  HISAT2( reads_ch)

}

