#! /usr/bin/env nextflow

picard_container = ''

process SamToFastq_help {
  label 'picard'

  container = "$picard_container"

  output: path 'picard-SamToFastq.txt'

  """
  picard SamToFastq > picard-SamToFastq.txt
  """
}

process SamToFastq_run {
    tag "$read_marked"
    label 'picard'
    publishDir "${params.outdir}/SamToFastq", mode: 'copy'

    input:
    path read_marked

    output:
    path "${read_marked.baseName}_samtofastq_interleaved.fq"

    script:
    """
    picard SamToFastq \
      I=${read_marked} \
      FASTQ=${read_marked.baseName}_samtofastq_interleaved.fq \
      CLIPPING_ATTRIBUTE=XT \
      CLIPPING_ACTION=2 \
      INTERLEAVE=true \
      NON_PF=true
    """
}
