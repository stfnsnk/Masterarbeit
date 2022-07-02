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
  
  container "${params.container_hisat}"
  publishDir "$params.out_dir/hisat/"

  input:
  path(reads)

  output:
  path '*.sam', emit: hisat2_sam_out_ch

// worked, can be optimised with sort?! samtools view -bS ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam > ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.bam
// the -p 8 flag didnt help
  script:
  """
  INDEX=`find -L ${params.hisat2_index}/ -name "*.1.ht2" | sed 's/.1.ht2//'` 
  hisat2 -x \$INDEX \
        -U $reads \
        -S ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam &> ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.summary.txt
  """
}

process sam_to_sorted_bam {
  tag"$sam_files.baseName"

  container "${params.container_samtools}"
  publishDir "$params.out_dir/bam_files/", mode: 'move'

  input:
  path(sam_files)

  output:
  path('*.bam')
  path('*.bai')

  script:
  """
  samtools flagstat ${sam_files} > ${sam_files.baseName}.flagstat.txt
  samtools sort -o ${sam_files.baseName}_sorted.bam ${sam_files}
  samtools index ${sam_files.baseName}_sorted.bam  
  """
}


workflow {

  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz}", checkIfExists: true)
    //.view()
    .set{reads_ch}

  HISAT2(reads_ch)
  sam_to_sorted_bam(HISAT2.out)
}

