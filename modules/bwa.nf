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
    tag "$name"
//    label 'process_medium'
    label 'bwa'
    publishDir "${params.outdir}/bwa", mode: 'copy'

//    when:
//    !params.skipQC && !params.skipFastQC

    input:
//    tuple val(name), file(reads) //from raw_reads_fastqc
    path(genome_fasta)

    output:
    path "$genome_fasta"
    path "$genome_fasta*"

    script:
    """
    bwa index $genome_fasta
    """
}
