#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "$baseDir/*.fastq.gz"
params.hisat2_index = "$baseDir/ref_genome_grcm38/grcm38_snp_tran.tar.gz"
//params.hisat2_indexfiles ="$baseDir/ref_genome_grcm38/grcm38_snptran"

log.info """\
         reads:        ${params.reads}
         hisat2_index: ${params.hisat2_index}
         """
         .stripIndent()

process UNCOMPRESS_GENOTYPE_INDEX {

  tag "$hisat2_index_ch.baseName"

  input:
  path(hisat2_index_ch)

  output:
  path('*')

  script:
  """
  tar xvzf $hisat2_index_ch
  """

}

process HISAT2 {

  tag "$reads.baseName"
  //container 'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3'

  input:
  path(hisat2_indices)
  path(reads)

  output:
  path('*.sam')

  script:
  """
  hisat2 -x $hisat2_indices \
         -U $reads \
         > ${reads.baseName}.sam
  """

}

workflow {

  Channel
    .fromPath( params.reads )
    .ifEmpty { exit 1, "reads was empty - no input files supplied" }
    .set { reads_ch }

//   Channel.fromPath( params.hisat2_index )
//     .ifEmpty { exit 1, "hisat2_index was empty - no input file supplied" }
//     .set { hisat2_index_ch }

  Channel
    .fromPath( params.hisat2_index )
    


 // UNCOMPRESS_GENOTYPE_INDEX(hisat2_index_ch)

  HISAT2(UNCOMPRESS_GENOTYPE_INDEX.out, reads_ch)

  //Command error:
  //(ERR): "grch38" does not exist
  //Exiting now ...

}

