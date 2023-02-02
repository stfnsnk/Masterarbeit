
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



//--------NOT IN USE--------//

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
