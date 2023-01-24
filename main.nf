#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.reads = "$baseDir/*.fastq.gz"
//params.hisat2_indexfiles ="$baseDir/ref_genome_grcm38/grcm38_snptran"



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

process featureCounts {
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
                -t exon \
                -g gene_id \
                --extraAttributes gene_name \
                -o ${bam_files.baseName}.counts.tsv ${bam_files} 
  """
}

process featureCounts_paired {
  
  // -p for paired end input
  //--extraAttributes gene_name  ------------noch nicht getestet

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
                -t exon \
                -g gene_id \
                --extraAttributes gene_name \
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
  // -k14 (lt. ONT Protokoll)
  //-Q	Ignore base quality in the input file 
  // -u finding canonical splice site b=both, f=transcript strand (lt. Paper "The long and the short of it; Dong X, Tian L et. al.")
  // --secondary=no TEST ob secondary alignment die featureCount anzahl verfälscht?
  // -p 1.0 (sollte bei --secondary=no eigentlich nicht nötig sein?!)
  tag "$sample_id"
  cpus 2 //for 6 core cpu with 64Gb RAM and queue = 2
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
           -k14 \
           -uf \
           --secondary=no \
           ${minimap2_index} \
           $input_files > ${sample_id}.sam   
  """
}

process pycoQC {//in Pipeline nicht einsetzbar < 64Gb RAM

  container "${params.container_pycoQC}"
  publishDir "${params.out_dir}/Nextflow_output/pycoQC/", mode: "move"
  
  input:
  path(input_seq_sum)       // guppy sequencing summary 
  path(input_barcode_sum)   // guppy_barcoder summary
  path(input_bam)           // bam files, optional if alignment metrics should additionally reported

  output:
  path("*")

  script:
  """
  if [ -e $input_bam ]; then 
    pycoQC  -f ${input_seq_sum} \
            -b ${input_barcode_sum} \
            -a ${input_bam}/*.bam \
            -o pycoQC_output.html \
            -j pycoQC_output.json
  else
    pycoQC  -f ${input_seq_sum} \
          -b ${input_barcode_sum} \
          -o pycoQC_output.html \
          -j pycoQC_output.json
  fi  
  """
}

process multiqc {

  container "$params.container_multiqc"
  containerOptions "--volume ${params.out_dir}/Nextflow_output:/multiqc_workdir"
  publishDir "${params.out_dir}/Nextflow_output/multiqc/", mode: "move"
  
  input:
  path(input_files)
  
 
  script:
  """
  cd /multiqc_workdir
  multiqc .
  """
}

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

//--------------------------------------------------------------------------------------//
//-----------------------------------WORKFLOW SECTION-----------------------------------//
//--------------------------------------------------------------------------------------//



workflow ont_pipeline {
  log.info """\
        =============================================
        ============DGE Analysis Pipeline============
        ============== ONT LONG-READ ================
        =============================================
        Author: Stefan Senk - BIOINF22
        reads: ${params.reads_folder}
        reference genome: ${params.minimap2_index_folder}
        features (GTF): ${params.features_to_count_folder}
        sample info file: ${params.sample_info}
        """
        .stripIndent()

  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
    .view()
    .set{reads_ch}

  Channel
    .fromPath( "${params.minimap2_index_folder}/*.fa*" )
    .collect()
    //.view()
    .set { minimap_ref_ch } 

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }


  minimap2(reads_ch, minimap_ref_ch)

  sam_to_sorted_filtered_bam(minimap2.out.minimap2_sam, minimap_ref_ch)

  //pycoQC("${params.pyco_seq_summary}", "${params.pyco_barcode_summary}", sam_to_sorted_filtered_bam.out.out_bamfiles.collect()) //noch nicht getestet!

  featureCounts_long(sam_to_sorted_filtered_bam.out.out_bamfiles, gene_annotation_file_ch)

  multiqc(featureCounts_long.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_long.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}

workflow short_read_pipeline_paired {
  log.info """\
         =============================================
         ============DGE Analysis Pipeline============
         ============PAIRED END SHORT-READ============
         =============================================
         Author: Stefan Senk - BIOINF22
         reads: ${params.reads_folder}
         reference genome: ${params.hisat2_index_folder}
         features (GTF): ${params.features_to_count_folder}
         sample info file: ${params.sample_info}
         """
         .stripIndent()

  Channel
    .fromFilePairs( "${params.reads_folder}/*_{R1,R2}*f*q.gz")
    .ifEmpty { exit 1, "no fastq files found at given path" }
    .set{ read_pair_ch }
  
  Channel
    .fromPath( "${params.hisat2_index_folder}/*.ht2" )
    .collect() 
    .map{ tuple( it.first().simpleName, it)}
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
 
  multiqc(featureCounts_paired.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_paired.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}

workflow short_read_pipeline_single {
  log.info """\
         =============================================
         ============DGE Analysis Pipeline============
         ============SINGLE END SHORT-READ============
         =============================================
         Author: Stefan Senk - BIOINF22
         reads: ${params.reads_folder}
         reference genome: ${params.hisat2_index_folder}
         features (GTF): ${params.features_to_count_folder}
         sample info file: ${params.sample_info}
         """
         .stripIndent()

  println("read files:")
  Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file)}
    .view()
    .set{reads_ch}
  
  Channel
    .fromPath( "${params.hisat2_index_folder}/*.ht2" )
    .collect() 
    .map{ tuple( it.first().simpleName, it)}
    .ifEmpty { exit 1, "no hisat2 index files found - path was empty" }
    .set { hisat2_index_ch }

  fastp_singleend(reads_ch)

  hisat2_singleend(fastp_singleend.out.trimmed_fastq, hisat2_index_ch)
  
  sam_to_sorted_bam(hisat2_singleend.out.sam_out_ch)

  featureCounts(sam_to_sorted_bam.out.sorted_bamfiles, "${params.features_to_count_folder}/*.gtf*")

  multiqc(featureCounts.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}



workflow filter_sam_and_count {
  Channel
    .fromPath("${params.input_path}/*.sam")
    .view()
    .set{ sam_input_ch }    

  Channel
    .fromPath( "${params.minimap2_index_folder}/*.fa*" )
    .collect()
    .view()
    .set { minimap_ref_ch } 

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }

  sam_to_sorted_filtered_bam(sam_input_ch, minimap_ref_ch)

  pycoQC("${params.pyco_seq_summary}", "${params.pyco_barcode_summary}", sam_to_sorted_filtered_bam.out.out_bamfiles.collect()) //noch nicht getestet!

  featureCounts_long(sam_to_sorted_filtered_bam.out.out_bamfiles, gene_annotation_file_ch)

  multiqc(featureCounts_long.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_long.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}

workflow DGE_analysis {
  Channel
    .fromPath("${params.DESEQ2_script}")
    .set{ r_script_ch }
    

   Channel
    .fromPath("${params.out_dir}/Nextflow_output/featureCounts", type: "dir" )
    .set{ counts_ch }
 

  Channel
    .fromPath("${params.input_path}")
    .set{ counts_ch }


  Channel
    .fromPath("${params.sample_info}")
    .set{ sample_info_ch }
    

  deseq2(sample_info_ch, counts_ch, r_script_ch)
}

workflow count_QC_DESEQ {

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }

  Channel
    .fromPath("${params.input_path}/*.bam")
    .set{ bam_to_count_ch }

  featureCounts_paired(bam_to_count_ch, gene_annotation_file_ch)
  
  multiqc(featureCounts_paired.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_paired.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}
