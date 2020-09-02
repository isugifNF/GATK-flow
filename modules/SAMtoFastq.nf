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
    tag "$name"
    label 'picard'
    publishDir "${params.outdir}/SamToFastq", mode: 'copy'

    input:
    val readname

    output:
    path "readname_samtofastq_interleaved.fq"

    script:
    """
    picard SamToFastq \
      I=${readname}_markilluminaadapters.bam \
      FASTQ=${readname}_samtofastq_interleaved.fq \
      CLIPPING_ATTRIBUTE=XT \
      CLIPPING_ACTION=2 \
      INTERLEAVE=true \
      NON_PF=true
    """
}
