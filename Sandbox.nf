#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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

process sam_to_sorted_filtered_bam {
  //-bh output in bam and with header
  // -T 
  /* -F 2324 exclude all reads with following flags
        read unmapped (0x4)
        read reverse strand (0x10)
        not primary alignment (0x100)
        supplementary alignment (0x800)
        2308 = reverse strand dabei

          samtools sort -o ${sam_files.baseName}_sorted.bam ${sam_files}
  samtools index ${sam_files.baseName}_sorted.bam
  */
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
                              -F 2324 | samtools sort -o ${sam_files.baseName}.sorted.filt.bam
  samtools index ${sam_files.baseName}.sorted.filt.bam 
  samtools flagstat ${sam_files.baseName}.sorted.filt.bam > ${sam_files.baseName}.sorted.filt.flagstat.txt  
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
  //  -Q	Ignore base quality in the input file 
  // -u finding canonical splice site b=both, f=transcript strand (lt. Paper "The long and the short of it; Dong X, Tian L et. al.")
  // --secondary=no TEST ob secondary alignment die featureCount anzahl verfälscht?

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


process minimap2_OLD {
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



process fastp_OLD {
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




//------------------------------------ONT TRIMMING TEST--------------------------------------//
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

workflow guppy_barcoder {

  //ONT Pipeline

   Channel 
    .fromPath("${params.reads_folder}/*.{fastq,fastq.gz,fq.gz,fq}", checkIfExists: true)
    .map { file -> tuple(file.simpleName, file)}
    .view()
    .set{reads_ch}



  guppy_barcode_trimming(reads_ch)

}
//-------------------------------------------------------------------------------------------//




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
  
  pycoQC_aln(pyco_input_seq_sum.collect(), pyco_input_bam, pyco_input_bai.collect())
}


//-----------------------------------Additional Workflows-----------------------------------//
//----------------------ADDITIONAL CONFIGURATION IS NEEDED BEFORE USE-----------------------//
//----------------------------params.input_path is your friend------------------------------//

workflow sam_to_DESEQ_long {
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

  //pycoQC("${params.pyco_seq_summary}", "${params.pyco_barcode_summary}", sam_to_sorted_filtered_bam.out.out_bamfiles.collect()) //noch nicht getestet!

  featureCounts_long(sam_to_sorted_filtered_bam.out.out_bamfiles, gene_annotation_file_ch)

  multiqc(featureCounts_long.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_long.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}

workflow filter_sam_to_DGE_short {//REV Filter makes no sense?!
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

  sam_to_sorted_filtered_bam_short(sam_input_ch)

  featureCounts_paired(sam_to_sorted_filtered_bam_short.out.out_bamfiles, gene_annotation_file_ch)

  //multiqc(featureCounts_paired.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_paired.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

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

workflow count_QC_DESEQ_paired {

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }

  Channel
    .fromPath("${params.input_path}/*.bam")
    .set{ bam_to_count_ch }

  featureCounts_paired(bam_to_count_ch, gene_annotation_file_ch)
  
  //multiqc(featureCounts_paired.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_paired.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}

workflow count_QC_DESEQ_long {

  Channel
    .fromPath( "${params.features_to_count_folder}/*.gtf*" )
    .collect()
    .set{ gene_annotation_file_ch }

  Channel
    .fromPath("${params.input_path}/*.bam")
    .set{ bam_to_count_ch }

  featureCounts_long(bam_to_count_ch, gene_annotation_file_ch)
  
  //multiqc(featureCounts_paired.out.featureCounts_summary_files.collect())

  deseq2("${params.sample_info}", featureCounts_long.out.count_files_ch.toList(), "${params.DESEQ2_script}") 

}


workflow NanoPlot { //getestet mit 1 und 3 bam files gleichzeitig - hängt nach einger Zeit, bricht jedoch nicht ab
  Channel
    .fromPath("${params.pyco_bam_files}")
    .view()
    .set{ pyco_input_bam }
  /*
  Channel
    .fromPath("${params.pyco_bam_files}/*.bai")
    .set{ pyco_input_bai }
  pyco_input_bai.collect()
  */
  NanoPlot(pyco_input_bam)
}

workflow {

 //download_mouse_reference()
 
 }