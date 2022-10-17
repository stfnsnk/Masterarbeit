#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process Trim_galore {
  tag "$input_files.baseName"

  container "$params.container_TrimGalore"
  publishDir "$params.out_dir/QC/"

  input:
  path(input_files)

  output:
  path("*.html")
 //-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
  script:
  """
  trim_galore  \
              --fastqc \
              $input_files
  """
  }



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
  featureCounts -a ${features} -t exon -g gene_id \
  -o ${bam_files.baseName}.counts.txt ${bam_files} 
  """
}

process minimap2 {
  tag "$input_files.baseName"

  container "$params.container_minimap2"
  publishDir "$params.out_dir/Nextflow_output/minimap2/"

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

process download_mouse_reference {

//TODO: hisat2 Index wird als komische datei runtergeladen, gunzip funktioniert nicht, muss manuell mit WinRAR?! entpackt werden

  script:
  """
  cd $projectDir
  wget -nc -O grcm38_snp_tran.gz -P reference_genome_GRCm38/Hisat_Index ${params.hisat2_index_url}
  wget -nc -P reference_genome_GRCm38/Gene_annotation ${params.GRCm38_gtf_url}
  wget -nc -P reference_genome_GRCm38/Minimap2_reference ${params.minimap2_GRCm38_ref_url} 
  """
}

process setup {

    script:
    """
    cd $launchDir
    mkdir -p LBIO_RNA_pipeline/INPUT/ONT_reads/Condition_A \
             LBIO_RNA_pipeline/INPUT/ONT_reads/Condition_B \
             LBIO_RNA_pipeline/INPUT/Illumina_reads/Condition_A \
             LBIO_RNA_pipeline/INPUT/Illumina_reads/Condition_B  
    """
}

workflow {


/*       Channel
        .fromPath("${params.reads_folder}/*.{fastq,fastq.gz}", type: 'file', checkIfExists: true)
        .view()
        .set{ reads_ch }
 */
//Trim_galore(reads_ch)

 /*     Channel
        .fromPath("${params.reads_folder}/*.sam", type: 'file', checkIfExists: true)
        .view()
        .set{ sam_files_ch }
 */
    //minimap2(reads_ch)

/*     Channel 
        .fromPath("${params.bam_files_to_count}/*.bam", type: 'file', checkIfExists: true)
        .view()
        .set{bam_file_ch}

    featureCounts(bam_file_ch, params.features_to_count)
 
 */
  download_mouse_reference()
 }