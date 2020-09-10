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
    tag "$sorted_ref"
    label 'samtools'
    publishDir "${params.outdir}/faidx", mode: 'copy'

    input:
    path sorted_ref

    output:
    path "${sorted_ref.baseName}.fasta.fai"

    script:
    """
    samtools faidx ${sorted_ref.baseName}.fasta
    """
}
