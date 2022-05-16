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
  publishDir "$params.out_dir/hisat/", mode: 'move'

  input:
  path(reads)

  output:
  path('*.sam')

  script:
  """
  hisat2 -x $params.hisat2_index \
         -U $reads \
         > ${reads.name.replaceAll('.fastq.gz','')}.sam
  """

}

workflow {

  Channel 
    .fromPath("${params.reads}/*.{fastq,fastq.gz}", checkIfExists: true)
    .set{reads_ch}
    
  HISAT2(reads_ch)

}

