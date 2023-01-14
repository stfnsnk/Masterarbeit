#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process trim_galore {
  tag "$sample_id"

  container "$params.container_TrimGalore"
  publishDir "$params.out_dir/trim_galore/"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*.html")
 
  script:
  """
  trim_galore --paired \
              --fastqc \
              $input_files
  """
}


process fastp_long_singleend {
  //following flags are set for ONT long reads
  // -Q, --disable_quality_filtering (quality filtering is enabled by default. If this option is specified, quality filtering is disabled)  
  // -L, --disable_length_filtering (length filtering is enabled by default. If this option is specified, length filtering is disabled)
  // --adapter_fasta (specify a FASTA file to trim reads by all the sequences in this FASTA file (string ([=auto])
  
  tag "$sample_id"

  container "${params.container_fastp}" 
  publishDir "${params.out_dir}/Nextflow_output/fastp/" , mode: "copy"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*.json")
  path("*.html")
  path("*_trimmed.fastq.gz") , emit: trimmed_fastq


  script:
  """
  fastp -i ${input_files[0]} \
        -o "${sample_id}_trimmed.fastq.gz" \
        --adapter_fasta
        --trim_poly_x \
        --disable_quality_filtering \
        --disable_length_filtering \
        --html "${sample_id}.html" \
        --json "${sample_id}.json" \
        --report_title "fastp report: $sample_id" \
        --overrepresentation_analysis
  """
}


process fastqc {
  tag "$sample_id"

  container "$params.container_TrimGalore"
  publishDir "$params.out_dir/FASTQC/"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*.html"), mode: move
  path("*.{fastq.gz,fastq,fq.gz,fq}"), mode: move

  script:
  """
  fastqc $input_files
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
  path('*.tsv')
  path('*.summary'), emit: featureCounts_summary_files

  script:
  """
  featureCounts -a ${features} \
                -t exon \
                -g gene_id \
                -o ${bam_files.baseName}.counts.tsv ${bam_files} 
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

  //TODO: TEST HISAT Index Download and tar
  script:
  """
  cd $projectDir
  mkdir -p reference_genome_GRCm38/Hisat_index
  [ -d $projectDir/reference_genome_GRCm38/Hisat_index/grcm38_snp_tran ] && wget -P reference_genome_GRCm38/Hisat_Index -O download.tar.gz ${params.hisat2_index_url} && tar -xzf download.tar.gz 
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

process fastp {
  tag "$sample_id"

  container "quay.io/biocontainers/fastp:0.23.2--h5f740d0_3" 
  publishDir "$params.out_dir/Nextflow_output/fastp/" , mode: "move"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*.html")
  path("*trimmed_fastq.gz") , emit: trimmed_fastq


  script:
  """
  fastp -i ${input_files[0]} \
        -o "${sample_id}_R1_trimmed.fastq.gz" \
        -I ${input_files[1]} \
        -O "${sample_id}_R2_trimmed.fastq.gz" \
        -g \
        -x \
        --html "$sample_id".html \
        --detect_adapter_for_pe
  """
}


process merge_fastq {
  tag "$sample_id"

  publishDir "$params.out_dir/fastq_merged" , mode: "move"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*")

  script:
  """
  cd $launchDir
  
  for read_nr in [1,2]; \
    do \
    for i in 'ls *_S?*_R${read_nr}*.fastq.gz'; \
      do  \
        head $i \
      done \
    done
  """

}

process pycoQC {
  cpus 6
  container "${params.container_pycoQC}"
  publishDir "${params.out_dir}/Nextflow_output/pycoQC/", mode: "move"
  
  input:
  path(input_seq_sum)       //guppy sequencing summary 
  path(input_barcode_sum)   // guppy_barcoder summary
  path(input_bam)           //directory containing bam files, optional if alignment metrics should additionally reported

  output:
  path("*")

  script:
  """
  if [ -e $input_bam ]; then 
    pycoQC  -f ${input_seq_sum} \
            -b ${input_barcode_sum} \
            -a ${input_bam}\
            -o pycoqc_output.html \
            -j pycoqc_output.json
  else
    pycoQC  -f ${input_seq_sum} \
          -b ${input_barcode_sum} \
          -o pycoqc_output.html \
          -j pycoqc_output.json
  fi  
  """
}

process pycoQC_aln {
  cpus 6
  container "${params.container_pycoQC}"
  publishDir "${params.out_dir}/Nextflow_output/pycoQC/", mode: "move"
  
  input:
  path(input_seq_sum)       //guppy sequencing summary 
  path(input_bam)           //directory containing bam files, optional if alignment metrics should additionally reported
  path(input_bai)

  output:
  path("*")

  script:
  """
 pycoQC -f ${input_seq_sum} \
        -a ${input_bam}\
        -o ${input_bam.simpleName}_pycoqc.html \
        -j ${input_bam.simpleName}_pycoqc.json
  """
}

process multiqc {

  container "${params.container_multiqc}"
  publishDir "${params.out_dir}/Nextflow_output/multiqc/", mode: "move"
  
  input:
  path(input_path)

  output:
  path("*")

  script:
  """
  cd ${input_path}
  multiqc .
  """
}

process guppy_barcode_trimming { 
  tag "$sample_id"

  container "${params.container_guppy_cpu}" 
  publishDir "${params.out_dir}/Nextflow_output/guppy/" , mode: "copy"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*")

  //--config dna_r9.4.1_450bps_sup_prom.cfg \
  script:
  """
 guppy_barcoder --input_path ${params.reads_folder} \
  --save_path \$PWD \
  --barcode_kits "EXP-NBD104 EXP-NBD114"
  """
}


workflow ont_pipeline {

  //ONT Pipeline

   Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file)}
    .view()
    .set{reads_ch}



  guppy_barcode_trimming(reads_ch)

}

workflow test {
  Channel
    .fromPath("${params.input_path}/*")
    .map { file -> tuple(file.simpleName, file)}
    .view()
    .set{reads_ch}


}

workflow featureCount_merge {
  Channel
    .fromPath("${params.input_path}")
    .collect()

}

workflow ONT_QC {
  Channel
    .fromPath("${params.pyco_seq_summary}")
    .collect()
    .set{ pyco_input_seq_sum }

  Channel
    .fromPath("${params.pyco_barcode_summary}")
    .collect()
    .set{ pyco_input_barcode_sum }
  
 
  pycoQC(pyco_input_seq_sum, pyco_input_barcode_sum)
}

workflow ONT_QC_aln {
  Channel
    .fromPath("${params.pyco_seq_summary}")
    .collect()
    .set{ pyco_input_seq_sum }
 
  Channel
    .fromPath("${params.pyco_bam_files}/*.bam")
    .view()
    .set{ pyco_input_bam }
  
  Channel
    .fromPath("${params.pyco_bam_files}/*.bai")
    .set{ pyco_input_bai }
  
  pycoQC_aln(pyco_input_seq_sum.collect(), pyco_input_bam, pyco_input_bai)
}

workflow {

/*  Channel
     .fromFilePairs("${params.reads_folder}/*_L00?_R{1,2}*.fastq*", checkIfExists: true)
    .view()
    .set{ read_pair_ch }
 */
  //fastp(read_pair_ch)

/*   Channel
  .fromFilePairs("${params.reads_folder}/*Ko1_{R1,R2}*fastq.gz", checkIfExists: true)
  .view()
  .set{ read_pair_ch }
 */

 //download_mouse_reference()
 
 //multiqc("${params.out_dir}/Nextflow_output/")

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



 }