#! /usr/bin/env nextflow

bioawk_container = ''

process sortSeq_help {
  label 'bioawk'

  output: path 'sortSeq_out.txt'

  """
  bioawk --help > sortSeq_out.txt
  """
}

process sortSeq_run {
    tag "${genome.fileName}"
    label 'bioawk'
    publishDir "${params.outdir}/sortSeq", mode: 'copy'

    input:
    val genome

    output:
    path "${genome.simpleName}_sorted.fasta"

    script:
    """
    bioawk -c fastx '{print}' $genome |\
     sort -k1,1V | awk '{print ">"\$1;print \$2}' |\
     fold > ${genome.simpleName}_sorted.fasta
    """
}

process seqLength_help {
  label 'bioawk'

  output: path 'seqLength_out.txt'

  """
  bioawk --help > seqLength_out.txt
  """
}

process seqLength_run {
    tag "$sorted_ref"
    label 'bioawk'
    publishDir "${params.outdir}/seqLength", mode: 'copy'

    input:
    path sorted_ref

    output:
    path "${sorted_ref.simpleName}.length"

    script:
    """
    bioawk -c fastx '{print \$name"\t"length(\$seq)}' $sorted_ref >\
     ${sorted_ref.simpleName}.length
    """
}
