
process deseq2 {
  
  tag "differential gene expression analysis with DESEQ2"
  
  container "$params.container_deseq2"
  
  publishDir "${params.out_dir}/Nextflow_output/DESEQ2/", mode: "move"


  input:
  path( sample_info_file )
  path( feature_count_files )
  path( script )

  output:
  path("*")

  script:
  """
  Rscript $script $sample_info_file $feature_count_files
  """
}
