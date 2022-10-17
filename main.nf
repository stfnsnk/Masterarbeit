#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.reads = "$baseDir/*.fastq.gz"
//params.hisat2_indexfiles ="$baseDir/ref_genome_grcm38/grcm38_snptran"

log.info """\
         reads: ${params.reads_folder}
         hisat2_index: ${params.hisat2_index}
         features (GTF): ${params.features_to_count}
         """
         .stripIndent()


process HISAT2 {

  tag "$reads.baseName"
  
  container "${params.container_hisat}"
  publishDir "$params.out_dir/Nextflow_output/hisat/"

  input:
  path(reads)

  output:
  path '*.sam' , emit: HISAT2_sam_out_ch
  path '*.txt'

  script:
  """
  INDEX=`find -L ${params.hisat2_index}/ -name "*.1.ht2" | sed 's/.1.ht2//'` 
  hisat2 --rna-strandness R \
          -x \$INDEX \
          -U $reads \
          -S ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam &> ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.summary.txt
  """
}

process sam_to_sorted_bam {
  tag"$sam_files.baseName"

  container "${params.container_samtools}"
  publishDir "$params.out_dir/Nextflow_output/bam_files/", mode: 'copy'

  input:
  path(sam_files)

  output:
  path('*.bam') , emit: sorted_bamfiles_ch
  path('*.bai')
  path('*.txt')

  script:
  """
  samtools flagstat ${sam_files} > ${sam_files.baseName}.flagstat.txt
  samtools sort -o ${sam_files.baseName}_sorted.bam ${sam_files}
  samtools index ${sam_files.baseName}_sorted.bam  
  """
}

process HISAT2_to_bam {

  tag "$reads.baseName"
  
  container "${params.container_HISAT2samtools}"
  publishDir "$params.out_dir/Nextflow_output/hisat/"

  input:
  path(reads)

  output:
  path '*.bam' , emit: HISAT2_bam_out_ch
  path '*.txt'

  script:
  """
  INDEX=`find -L ${params.hisat2_index}/ -name "*.1.ht2" | sed 's/.1.ht2//'` 
  hisat2 --rna-strandness R \
        -x \$INDEX \
        -U $reads \
        -S ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam &> ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.hisat_summary.txt
  samtools flagstat ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam > ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.flagstat.txt
  samtools sort -o ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}_sorted.bam ${reads.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam
  """
}

process featureCounts {
  tag "$bam_files.baseName"

  container "$params.container_subread"
  publishDir "$params.out_dir/Nextflow_output/featureCounts/", mode: "move"

  input:
  path(bam_files)
  path(features)

  output:
  path('*.txt')
  path('*.summary')

  script:
  """
  featureCounts -a ${features} \
                -s 2 \
                -t exon \
                -g gene_id \
                -o ${bam_files.baseName}.counts.txt ${bam_files} 
  """
}

process minimap2 {
  tag "$input_files.baseName"

  container "$params.container_minimap2"

  input:
  path(input_files)

  output:
  path("*.sam"), emit: minimap2_sam

  script:
  """
  minimap2 -ax splice \
            ${params.minimap2_index} \
            $input_files > ${input_files.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam   
  """
}
workflow {

  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    //.view()
    .set{reads_ch}

  HISAT2(reads_ch)
  sam_to_sorted_bam(HISAT2.out.HISAT2_sam_out_ch)

  //minimap2(reads_ch)
  //sam_to_sorted_bam(minimap2.out.minimap2_sam)
  
  //HISAT2_to_bam(reads_ch)
  //featureCounts(sam_to_sorted_bam.out.sorted_bamfiles_ch, params.features_to_count)
}

