#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process hisat_version {
    
    container "${params.container_hisat}"

    script:
    """
    hisat2 --version
    """

}

process test {

    tag "${reads.name.replaceAll('.fastq.gz','')}"
    publishDir "$params.out_dir/test/", mode: 'move'

    input:
    path reads

    output:
    path('*.txt')
    
    script:
    """
    zcat $reads | head -1 > firstlines_${reads.name.replaceAll('.fastq.gz','')}.txt
    """

}



workflow {

    println ("Hy starting my first workflow with" + params.container_hisat)
    println ("with the reads from " + params.reads)
    Channel
        .fromPath("${params.reads}/*.{fastq,fastq.gz}", checkIfExists: true)
        .view()
        .set{ reads_ch }

    test(reads_ch)
    //hisat_version()
    
}