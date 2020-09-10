#! /usr/bin/env nextflow

bwa_container = 'kathrinklee/bwa'

process bwa_help {
  label 'bwa'

  container = "$bwa_container"

  output: path 'bwa_out.txt'

  """
  echo "bwa" > bwa_out.txt
  bwa &>> bwa_out.txt
  """
}

process bwa_index {
    tag "$sorted_ref"
//    label 'process_medium'
    label 'bwa'
    publishDir "${params.outdir}/bwa", mode: 'copy'

//    when:
//    !params.skipQC && !params.skipFastQC

    input:
//    tuple val(name), file(reads) //from raw_reads_fastqc
    path(sorted_ref)

    output:
    path "$sorted_ref"
    path "$sorted_ref*"

    script:
    """
    bwa index $sorted_ref
    """
}
