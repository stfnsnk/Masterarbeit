#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.reads = "$baseDir/*.fastq.gz"
//params.hisat2_indexfiles ="$baseDir/ref_genome_grcm38/grcm38_snptran"

log.info """\
         reads: ${params.reads_folder}
         hisat2_index: ${params.hisat2_index}
         """
         .stripIndent()


process HISAT2 {

  tag "$reads.baseName"
  
  container "${params.container_HISAT2samtools}"
  publishDir "$params.out_dir/hisat/", mode: 'move'

  input:
  path(reads)

  output:
  path('*.sam')

  script:
  """
  samtools
  """


  // """
  // INDEX=`find -L ${params.hisat2_index}/ -name "*.1.ht2" | sed 's/.1.ht2//'` 
  
  // hisat2 -x \$INDEX \
  //        -U $reads \
  //        > ${reads.name.replaceAll('.fastq.gz','')}.sam

  // """
// samtools view ${reads.name.replaceAll('.fastq.gz','')}.sam > ${reads.name.replaceAll('.fastq.gz','')}.bam
// samtools flagstat ${reads.name.replaceAll('.fastq.gz','')}.sam > ${reads.name.replaceAll('.fastq.gz','')}.flagstat.txt


process samtools {

  tag "$reads.baseName"

  con

}


}

workflow {

  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz}", checkIfExists: true)
    //.view()
    .set{reads_ch}

  HISAT2(reads_ch)

}

