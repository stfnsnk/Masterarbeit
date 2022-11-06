#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.reads = "$baseDir/*.fastq.gz"
//params.hisat2_indexfiles ="$baseDir/ref_genome_grcm38/grcm38_snptran"

log.info """\
         reads: ${params.reads_folder}
         hisat2_index: ${params.hisat2_index}
         features (GTF): ${params.features_to_count_folder}
         """
         .stripIndent()


process fastp {
  tag "$sample_id"

  container "${params.container_fastp}" 
  publishDir "${params.out_dir}/Nextflow_output/fastp/" , mode: "copy"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*.html")
  tuple val("$sample_id"), path("*{R1,R2}_trimmed.fastq.gz") , emit: trimmed_fastq_pair


  script:
  """
  fastp -i ${input_files[0]} \
        -o "${sample_id}_R1_trimmed.fastq.gz" \
        -I ${input_files[1]} \
        -O "${sample_id}_R2_trimmed.fastq.gz" \
        --trim_poly_x \
        --html "$sample_id".html \
        --detect_adapter_for_pe \
        --overrepresentation_analysis
  """
}



process hisat2_singleend {

  tag "$reads.baseName"
  
  container "${params.container_hisat}"
  publishDir "$params.out_dir/Nextflow_output/hisat/"

  input:
  path(reads)

  output:
  path '*.sam' , emit: sam_out_ch
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


process hisat2_paired {

  tag "$sample_id"
  
  container "${params.container_hisat}"
  publishDir "${params.out_dir}/Nextflow_output/hisat/"

  input:
  tuple val(sample_id), path(input_files)
  path( index_files )

  output:
  path('*.sam') , emit: hisat2_sam_out_ch
  path('*.txt')

  script:
  """
  INDEX=`find -L ${index_files}/ -name "*.?.ht2" | sed 's/.?.ht2//'` 
  hisat2 --rna-strandness RF \
          -x \$INDEX \
          -1 ${input_files[0]} \
          -2 ${input_files[1]} \
          -S ${sample_id}.sam &> ${sample_id}.hisat_summary.txt
  """
}

process sam_to_sorted_bam {
  tag"$sam_files.baseName"

  container "${params.container_samtools}"
  publishDir "$params.out_dir/Nextflow_output/bam_files/", mode: 'copy'

  input:
  path(sam_files)

  output:
  path('*.bam') , emit: sorted_bamfiles
  path('*.bai')
  path('*.txt')

  script:
  """
  samtools flagstat ${sam_files} > ${sam_files.baseName}.flagstat.txt
  samtools sort -o ${sam_files.baseName}_sorted.bam ${sam_files}
  samtools index ${sam_files.baseName}_sorted.bam  
  """
}

process hisat2_to_bam {

  tag "$reads.baseName"
  
  container "${params.container_HISAT2samtools}"
  publishDir "$params.out_dir/Nextflow_output/hisat/"

  input:
  path(reads)

  output:
  path '*.bam' , emit: hisat2_bam_out_ch
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

process featureCounts_paired {
  
  // -p for paired end input

  tag "$bam_files.baseName"

  container "${params.container_subread}"
  publishDir "${params.out_dir}/Nextflow_output/featureCounts/", mode: "move"

  input:
  path(bam_files)
  path(features)

  output:
  path('*.tsv')
  path('*.summary')

  script:
  """
  featureCounts -a ${features} \
                -s 2 \
                -p \
                -t exon \
                -g gene_id \
                -o ${bam_files.baseName}.counts.tsv ${bam_files} 
  """
}

process featureCounts_long {
  
  //the -L flag turns on long-read counting mode

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
                -L \
                -t exon \
                -g gene_id \
                -o ${bam_files.baseName}.counts.txt ${bam_files} 
  """
}


process minimap2 {
  tag "$input_files.baseName"

  container "$params.container_minimap2"

  input:
  path( input_files )
  path( minimap_index )

  output:
  path("*.sam"), emit: minimap2_sam

  script:
  """
  minimap2 -ax splice \
            ${minimap2_index} \
            $input_files > ${input_files.name.replaceAll("['.fastq.gz'|'.fastq']",'')}.sam   
  """
}



workflow {

/*   Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    //.view()
    .set{reads_ch}
 */

   Channel
    .fromPath( params.hisat2_index )
    .collect() //
    .ifEmpty { exit 1, "no hisat2 index files found - path was empty" }
    .set { hisat2_index_ch }
 
  Channel
    .fromFilePairs( "${params.reads_folder}/*_{R1,R2}*f*q.gz")
    .ifEmpty { exit 1, "no fastq files found at given path" }
    .view()
    .set{ read_pair_ch }

  Channel
    .fromPath( params.minimap2_index )
    .collect()
    .set { minimap_index_ch }

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file }
 
  fastp(read_pair_ch)

  hisat2_paired(fastp.out.trimmed_fastq_pair, hisat2_index_ch)
  //hisat2_paired( read_pair_ch, hisat2_index_ch )

  sam_to_sorted_bam(hisat2_paired.out.sam_out_ch)

  //minimap2(reads_ch, minimap_index_ch)
  //sam_to_sorted_bam(minimap2.out.minimap2_sam)
  
  //featureCounts_paired(sam_to_sorted_bam.out.sorted_bamfiles, gene_annotation_file)
}

