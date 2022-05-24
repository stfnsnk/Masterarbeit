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

    tag "${reads.name.replaceAll('.fastq.gz','')}"
    publishDir "$params.out_dir/firstlines/", mode: 'move'

    input:
    path reads

    output:
    path('*.txt')
    
    script:
    """
    zcat $reads | head -1 > firstline_${reads.name.replaceAll('.fastq.gz','')}.txt
    """

}



workflow {

    println ("Hy starting my first workflow with " + params.container_hisat)
    println ("with the reads from " + params.reads_folder)
    Channel
        .fromPath("${params.reads_folder}/*.{fastq,fastq.gz}", type: 'file', checkIfExists: true)
        .view()
        .set{ reads_ch }

    first_line(reads_ch)
    //hisat_version()
    
}