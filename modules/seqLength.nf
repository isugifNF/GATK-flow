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
    tag "$name"
    label 'bioawk'
    publishDir "${params.outdir}/seqLength", mode: 'copy'

    input:
    val shortname
    
    output:
    path "shortname.length"

    script:
    """
    bioawk -c fastx '{print \$name"\t"length(\$seq)}' $genome > ${shortname}.length
    """
}
