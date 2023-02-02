
process fastp_singleend {
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
        --trim_poly_x \
        --html "${sample_id}.fastp.html" \
        --json "${sample_id}.fastp.json" \
        --report_title "fastp report: $sample_id" \
        --overrepresentation_analysis
  """
}

process fastp_paired {
  tag "$sample_id"

  container "${params.container_fastp}" 
  publishDir "${params.out_dir}/Nextflow_output/fastp/" , mode: "copy"

  input:
  tuple val(sample_id), path(input_files)

  output:
  path("*.json")
  path("*.html")
  tuple val("$sample_id"), path("*{R1,R2}_trimmed.fastq.gz") , emit: trimmed_fastq_pair


  script:
  """
  fastp -i ${input_files[0]} \
        -o "${sample_id}_R1_trimmed.fastq.gz" \
        -I ${input_files[1]} \
        -O "${sample_id}_R2_trimmed.fastq.gz" \
        --trim_poly_x \
        --html "${sample_id}.fastp.html" \
        --json "${sample_id}.fastp.json" \
        --report_title "fastp report: $sample_id" \
        --detect_adapter_for_pe \
        --overrepresentation_analysis
  """
}


//--------NOT IN USE--------//

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

