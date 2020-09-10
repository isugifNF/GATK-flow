#! /usr/bin/env nextflow

picard_container = ''

process createSeqDict_help {
  label 'picard'

  container = "$picard_container"

  output: path 'picard-seqCreateDict.txt'

  """
  picard CreateSequenceDictionary > picard-seqCreateDict.txt
  """
}

process createSeqDict_run {
    tag "$sorted_ref"
    label 'picard'
    publishDir "${params.outdir}/createSeqDict", mode: 'copy'

    input:
    path sorted_ref

    output:
    path "${sorted_ref.baseName}.dict"

    script:
    """
    picard CreateSequenceDictionary \
      REFERENCE=${sorted_ref} \
      OUTPUT=${sorted_ref.baseName}.dict
    """
}
