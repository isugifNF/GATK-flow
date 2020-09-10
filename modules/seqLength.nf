#! /usr/bin/env nextflow

bioawk_container = 'kathrinklee/bwa'

process seqLength_help {
  label 'bioawk'

  container = "$bioawk_container"

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
    path "${sorted_ref.baseName}.length"

    script:
    """
    bioawk -c fastx '{print \$name"\t"length(\$seq)}' $sorted_ref >\
     ${sorted_ref.baseName}.length
    """
}
