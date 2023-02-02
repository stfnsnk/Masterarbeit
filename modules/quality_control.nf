

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

process pycoQC {//may not work with big bam files and host computer with less than 64Gb RAM

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


//--------NOT IN USE--------//

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

process pycoQC_sum_only {
  cpus 6
  container "${params.container_pycoQC}"
  publishDir "${params.out_dir}/Nextflow_output/pycoQC/", mode: "move"
  
  input:
  path(input_seq_sum)       //guppy sequencing summary 
  path(input_barcode_sum)   // guppy_barcoder summary

  output:
  path("*")

  script:
  """
    pycoQC  -f ${input_seq_sum} \
            -b ${input_barcode_sum} \
            -o pycoqc_output.html \
            -j pycoqc_output.json
  
  """
}

process pycoQC_aln {
  container "${params.container_pycoQC}"
  publishDir "${params.out_dir}/Nextflow_output/pycoQC/", mode: "move"
  
  input:
  path(input_seq_sum)       //guppy sequencing summary to corresponding bam file
  path(input_bam)           //single bam file
  path(input_bai)

  output:
  path("*")

  script:
  """
 pycoQC -f ${input_seq_sum} \
        -a ${input_bam} \
        -o ${input_bam.simpleName}_pycoqc.html \
        -j ${input_bam.simpleName}_pycoqc.json
  """
}

process NanoPlot {//may be better suited for alignment QC integrated in the pipeline; does not need big sequencing_summary file - TEST failed!!
  container "${params.container_nanoplot}"
  publishDir "${params.out_dir}/Nextflow_output/NanoPlot/", mode: "move"
  
  input:
  path(input_bam)          //single bam file or List of bam files 
  path(input_bai)

  output:
  path("*")

  script:
  """
  NanoPlot  -t 12 \
            --color yellow \
            --bam $input_bam \
            --downsample 10000 \
            -o bamplots_downsampled
  """
}


