#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process hisat_version {
    
    container "${params.container_hisat}"

    script:
    """
    hisat2 --version
    """

}

process first_line {

    tag "${reads.name.replaceAll("['.fastq.gz'|'.fastq']","")}"
    publishDir "$params.out_dir/firstlines/", mode: 'move'

    input:
    path reads

    output:
    path('*.txt')
    
    script:
    """
    zcat $reads | head -1 > firstline_${reads.name.replaceAll("['.fastq.gz'|'.fastq']","")}.txt
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
  samtools sort -o ${sam_files.baseName}.bam ${sam_files}
  samtools index ${sam_files.baseName}.bam  
  samtools flagstat ${sam_files.baseName}.bam > ${sam_files.baseName}.flagstat.txt
  """
}



workflow {

    println ("Hy starting my first workflow with " + params.container_hisat)
    println ("with the reads from " + params.reads_folder)
    Channel
        .fromPath("${params.reads_folder}/*.{fastq,fastq.gz}", type: 'file', checkIfExists: true)
        .view()
        .set{ reads_ch }

    Channel
        .fromPath("${params.reads_folder}/*.sam", type: 'file', checkIfExists: true)
        .view()
        .set{ sam_files_ch }

    sam_to_sorted_bam(sam_files_ch)
    //first_line(reads_ch)
    //hisat_version()
    
}