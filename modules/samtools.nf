
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
  samtools sort -o ${sam_files.baseName}.sorted.bam ${sam_files}
  samtools flagstat ${sam_files.baseName}.sorted.bam > ${sam_files.baseName}.flagstat.txt
  samtools index ${sam_files.baseName}.sorted.bam  
  """
}

process sam_to_sorted_filtered_bam {
  //-bh output in bam and with header
  // -T 
  /* -F 2324 !exclude! all reads with following flags
        read unmapped (0x4)
        read reverse strand (0x10) => filter excludes all genes on reverse strand?!
        not primary alignment (0x100)
        supplementary alignment (0x800)
        2308 = reverse strand dabei
  */
  // flagstat of the sorted.filt.bam resulted in only primary mapping reads- TEST OK

  tag"$sam_files.baseName"

  container "${params.container_samtools}"
  publishDir "$params.out_dir/Nextflow_output/bam_files/", mode: 'copy'

  input:
  path( sam_files )
  path( ref_genome_fasta )

  output:
  path('*.sorted.filt.bam') , emit: out_bamfiles
  path('*.bai')
  path('*.txt')

  script:
  """
  samtools flagstat ${sam_files} > ${sam_files.baseName}.flagstat.txt
  samtools view ${sam_files}  -bh \
                              -t ${ref_genome_fasta} \
                              -F 2308 | samtools sort -o ${sam_files.baseName}.sorted.filt.bam
  samtools index ${sam_files.baseName}.sorted.filt.bam 
  """
}

process sam_to_sorted_filtered_bam_short {
  //-bh output in bam and with header
  // -T 
  /* -F 2324 !exclude! all reads with following flags
        read unmapped (0x4)
        read reverse strand (0x10)
        not primary alignment (0x100)
        supplementary alignment (0x800)
        2308 = reverse strand dabei
  */
  // flagstat of the sorted.filt.bam resulted in only primary mapping reads- TEST OK

  tag"$sam_files.baseName"

  container "${params.container_samtools}"
  publishDir "$params.out_dir/Nextflow_output/bam_files/", mode: 'copy'

  input:
  path( sam_files )

  output:
  path('*.sorted.filt.bam') , emit: out_bamfiles
  path('*.bai')
  path('*.txt')

  script:
  """
  samtools flagstat ${sam_files} > ${sam_files.baseName}.flagstat.txt
  samtools view ${sam_files}  -bh \
                              -F 2308 | samtools sort -o ${sam_files.baseName}.sorted.filt.bam
  samtools index ${sam_files.baseName}.sorted.filt.bam 
  """
}
