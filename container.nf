#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process hisat_version {
    
    container "${params.container_hisat}"

    script:
    """
    hisat2 --version
    """

}

workflow {

    println("Hy starting my first workflow with $params.container_hisat")

    hisat_version()
    print hisat_version.out
}