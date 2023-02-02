
process featureCounts {
  // -s < intorstring > (isStrandSpecific) => 2 means reverse; for PE reverse corresponds to orientation of first fragment/read

  tag "$bam_files.baseName"

  container "$params.container_subread"
  publishDir "$params.out_dir/Nextflow_output/featureCounts/", mode: "move"

  input:
  path(bam_files)
  path(features)

  output:
  path('*.tsv'), emit: count_files_ch
  path('*.summary'), emit: featureCounts_summary_files

  script:
  """
  featureCounts -a ${features} \
                -s 2 \
                -t exon \
                -g gene_id \
                --extraAttributes gene_name \
                -o ${bam_files.baseName}.counts.tsv ${bam_files} 
  """
}

process featureCounts_long {
  
  //the -L flag turns on long-read counting mode
  /* --primary (primaryOnly) If specified, only primary alignments will be counted. 
                          Primary and secondary alignments are identified using 
                          bit 0x100 in the Flag field of SAM/BAM files. 
                          All primary alignments in a dataset will be counted 
                          no matter they are from multimapping reads or not 
                          (ie. ‘-M’ is ignored). */
  // -O (AllowMultiOverlap) not recommended for RNA-Seq
  //−−extraAttributes gene_name  ------------noch nicht getestet
  
  tag "$bam_files.baseName"

  container "$params.container_subread"
  publishDir "$params.out_dir/Nextflow_output/featureCounts/", mode: "copy"

  input:
  path(bam_files)
  path(features) 

  output:
  path('*.tsv'), emit: count_files_ch
  path('*.summary'), emit: featureCounts_summary_files

  script:
  """
  featureCounts -a ${features} \
                -L \
                -t exon \
                -g gene_id \
                --primary \
                --extraAttributes gene_name \
                -o ${bam_files.simpleName}.counts.tsv ${bam_files} 
  """
}

process featureCounts_paired {
  // -s < intorstring > (isStrandSpecific) => 2 means reverse; for PE reverse corresponds to orientation of first fragment/read
  // -p for paired end input
  // --countReadPairs - should be specified according to manual if -p is set (since Release 2.0.2)

  tag "$bam_files.baseName"

  container "${params.container_subread}"
  publishDir "${params.out_dir}/Nextflow_output/featureCounts/", mode: "copy"

  input:
  path(bam_files)
  path(features)

  output:
  path('*.tsv'), emit: count_files_ch
  path('*.summary'), emit: featureCounts_summary_files

  script:
  """
  featureCounts -a ${features} \
                -s 2 \
                -p \
                --countReadPairs \
                -t exon \
                -g gene_id \
                --extraAttributes gene_name \
                -o ${bam_files.simpleName}.counts.tsv ${bam_files} 
  """
}
