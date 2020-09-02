#! /usr/bin/env nextflow

bioawk_container = ''

process sortSeq_help {
  label 'bioawk'

  container = "$bioawk_container"
  
  output: path 'sortSeq_out.txt'

  """
  bioawk --help > sortSeq_out.txt
  """
}

process sortSeq_run {
    tag "$name"
    label 'bioawk'
    publishDir "${params.outdir}/sortSeq", mode: 'copy'

    input:
    val genome
    val shortname
    
    output:
    path "shortname.fasta"

    script:
    """
    bioawk -c fastx '{print}' $genome | sort -k1,1V | awk '{print ">"\$1;print \$2}' | fold > ${shortname}.fasta 
    """
}
