#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.reads = "$baseDir/*.fastq.gz"
//params.hisat2_indexfiles ="$baseDir/ref_genome_grcm38/grcm38_snptran"

log.info """\
         reads: ${params.reads_folder}
         hisat2_index: ${params.hisat2_index_folder}
         features (GTF): ${params.features_to_count_folder}
         """
         .stripIndent()


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
        --html "${sample_id}.html" \
        --json "${sample_id}.json" \
        --report_title "fastp report: $sample_id" \
        --detect_adapter_for_pe \
        --overrepresentation_analysis
  """
}

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
        --html "${sample_id}.html" \
        --json "${sample_id}.json" \
        --report_title "fastp report: $sample_id" \
        --overrepresentation_analysis
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



process hisat2_singleend {

  tag "$reads.baseName"
  
  container "${params.container_hisat}"
  publishDir "${params.out_dir}/Nextflow_output/hisat/", pattern: "*.txt", mode: "move"

  input:
  tuple val(sample_id), path(input_files)
  tuple val( hisat_index_filename ), path( index_file_path )

  output:
  path( '*.sam' ) , emit: sam_out_ch
  path( '*.txt' )

  script:
  """
  hisat2 --rna-strandness R \
          --new-summary \
          --summary-file ${sample_id}.hisat_summary.txt \
          -x ${params.hisat2_index_folder}/${hisat_index_filename} \
          -U ${input_files} \
          -S ${sample_id}.sam
  """
}


process hisat2_paired {

  tag "$sample_id"
  
  container "${params.container_hisat}"
  publishDir "${params.out_dir}/Nextflow_output/hisat/", pattern: "*.txt", mode: "move"

  input:
  tuple val(sample_id), path(input_files)
  tuple val( hisat_index_filename ), path( index_file_path )

  output:
  path('*.sam') , emit: sam_out_ch
  path('*.txt')

  script:
  """
  hisat2 --rna-strandness RF \
          --new-summary \
          --summary-file ${sample_id}.hisat_summary.txt \
          -x ${params.hisat2_index_folder}/${hisat_index_filename} \
          -1 ${input_files[0]} \
          -2 ${input_files[1]} \
          -S ${sample_id}.sam 
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

process sam_to_sorted_filtered_bam {
  //-bh output in bam and with header
  // -T 
  /* -F 2324 exclude all reads with following flags
        read unmapped (0x4)
        read reverse strand (0x10)
        not primary alignment (0x100)
        supplementary alignment (0x800)
  */
  tag"$sam_files.baseName"

  container "${params.container_samtools}"
  publishDir "$params.out_dir/Nextflow_output/bam_files/", mode: 'copy'

  input:
  path(sam_files)
  path( ref_genome_fasta )

  output:
  path('*.bam') , emit: sorted_bamfiles
  path('*.bai')
  path('*.txt')

  script:
  """
  samtools flagstat ${sam_files} > ${sam_files.baseName}.flagstat.txt
  samtools view ${sam_files} -bh \
                                      -T ${ref_genome_fasta} \
                                      -F 2324 > ${sam_files.baseName}.filt.bam
  samtools sort -o ${sam_files.baseName}_filt_sorted.bam ${sam_files.baseName}.filt.bam
  samtools index ${sam_files.baseName}_filt_sorted.bam 
  samtools flagstat ${sam_files.baseName}_filt_sorted.bam > ${sam_files.baseName}_filt_sorted.flagstat.txt  
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
  path('*.summary'), emit: featureCounts_summary_files

  script:
  """
  featureCounts -a ${features} \
                -s 2 \
                -p \
                -t exon \
                -g gene_id \
                -o ${bam_files.simpleName}.counts.tsv ${bam_files} 
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
                -L \
                -t exon \
                -g gene_id \
                --primary \
                -o ${bam_files.simpleName}.counts.tsv ${bam_files} 
  """
}


process minimap2 {
  
  //-a output in SAM format
  // -x setting preset = 
  /* splice =Long-read spliced alignment
      (-k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -b0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes). 
      In the splice mode, 
      1) long deletions are taken as introns and represented as the ‘N’ CIGAR operator; 
      2) long insertions are disabled; 
      3) deletion and insertion gap costs are different during chaining; 
      4) the computation of the ‘ms’ tag ignores introns to demote hits to pseudogenes. */
  //-Q	Ignore base quality in the input file 
  // -u finding canonical splice site f=transcript strand (lt. Paper "The long and the short of it; Dong X, Tian L et. al.")
  // --secondary=no TEST ob secondary alignment die featureCount anzahl verfälscht?
  // -p 1.0 (sollte bei --secondary=no eigentlich nicht nötig sein?!)
  tag "$sample_id"
  cpus 3 //for 6 core cpu with 64Gb RAM and queue = 2
  container "$params.container_minimap2"

  input:
  tuple val(sample_id), path(input_files)
  path( minimap2_index )

  output:
  path("*.sam"), emit: minimap2_sam

  script:
  """
  minimap2 -a \
           -x splice \
           -uf \
           -p 1.0 \
           --secondary=no \
           ${minimap2_index} \
           $input_files > ${sample_id}.sam   
  """
}

process multiqc {

  container "$params.container_multiqc"
  containerOptions "--volume ${params.out_dir}/Nextflow_output:/multiqc_workdir"
  publishDir "${params.out_dir}/Nextflow_output/multiqc/", mode: "move"
  
  input:
  path(input_files).collect()
  

  output:
  path("*")

  script:
  """
  cd /multiqc_workdir \
  multiqc .
  """
}


workflow ont_pipeline {
  
  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file)}
    .view()
    .set{reads_ch}

  Channel
    .fromPath( "${params.minimap2_index_folder}/*.fa*" )
    .collect()
    .view()
    .set { minimap_index_ch } 

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }


  minimap2(reads_ch, minimap_index_ch)
  sam_to_sorted_bam(minimap2.out.minimap2_sam)
  featureCounts_long(sam_to_sorted_bam.out.sorted_bamfiles, gene_annotation_file_ch)

  //multiqc(featureCounts_long.out.featureCounts_summary_files)
}

workflow short_read_pipeline_paired {

  Channel
    .fromFilePairs( "${params.reads_folder}/*_{R1,R2}*f*q.gz")
    .ifEmpty { exit 1, "no fastq files found at given path" }
    .view()
    .set{ read_pair_ch }
  
  Channel
    .fromPath( "${params.hisat2_index_folder}/*.ht2" )
    .collect() 
    .map{ tuple( it.first().simpleName, it)}
    .view()
    .ifEmpty { exit 1, "no hisat2 index files found - path was empty" }
    .set { hisat2_index_ch }  

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }
   
  fastp_paired(read_pair_ch)

  hisat2_paired(fastp_paired.out.trimmed_fastq_pair, hisat2_index_ch)
  sam_to_sorted_bam(hisat2_paired.out.sam_out_ch)

  featureCounts_paired(sam_to_sorted_bam.out.sorted_bamfiles, gene_annotation_file_ch)
 
  //multiqc(featureCounts_paired.out.featureCounts_summary_files)
}

workflow short_read_pipeline_single {

  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file)}
    .view()
    .set{reads_ch}
  
  Channel
    .fromPath( "${params.hisat2_index_folder}/*.ht2" )
    .collect() 
    .map{ tuple( it.first().simpleName, it)}
    .view()
    .ifEmpty { exit 1, "no hisat2 index files found - path was empty" }
    .set { hisat2_index_ch }

fastp_singleend(reads_ch)

hisat2_singleend(fastp.out.trimmed_fastq, hisat_index_ch)
sam_to_sorted_bam(hisat2_singleend.out.sam_out_ch)

featureCounts(sam_to_sorted_bam.out.sorted_bamfiles, gene_annotation_file_ch)

//multiqc(featureCounts_paired.out.featureCounts_summary_files)

}
