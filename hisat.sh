#!/usr/bin/bash

docker run -v $PWD:/mnt/data quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3 hisat2 \
                -x /mnt/data/genome_snp_tran/genome_snp_tran \
                -U /mnt/data/D1_Ko_S4_L001_R1_001.fastq.gz \
                -S /mnt/data/eg1.sam