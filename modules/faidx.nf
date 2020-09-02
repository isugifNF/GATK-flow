#! /usr/bin/env nextflow

samtools_container = ''

process faidx_help {
  label 'samtools'

  container = "$samtools_container"
  
  output: path 'samtools-faidx.txt'

  """
  samtools faidx > samtools-faidx.txt
  """
}

process faidx_run {
    tag "$name"
    label 'samtools'
    publishDir "${params.outdir}/faidx", mode: 'copy'

    input:
    val shortname
    
    output:
    path "shortname.fasta.fai"

    script:
    """
    samtools faidx ${shortname}.fasta
    """
}
