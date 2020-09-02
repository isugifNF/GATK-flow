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
    tag "$name"
    label 'picard'
    publishDir "${params.outdir}/createSeqDict", mode: 'copy'

    input:
    val shortname
    
    output:
    path "shortname.dict"

    script:
    """
    picard CreateSequenceDictionary REFERENCE=${shortname}.fasta OUTPUT=${shortname}.dict
    """
}
