#! /usr/bin/env nextflow

picard_container = ''

process markAdapters_help {
  label 'picard'

  container = "$picard_container"
  
  output: path 'picard-markAdapters.txt'

  """
  picard MarkIlluminaAdapters > picard-markAdapters.txt
  """
}

process markAdapters_run {
    tag "$name"
    label 'picard'
    publishDir "${params.outdir}/markAdapters", mode: 'copy'

    input:
    val readname
    val TMPDIR
    
    output:
    path "readname_markilluminaadapters.bam"
    path "readname_markilluminaadapters_metrics.txt"

    script:
    """
    picard MarkIlluminaAdapters \
      I=${readname}_fastqtosam.bam \
      O=${readname}_markilluminaadapters.bam \
      M=${readname}_markilluminaadapters_metrics.txt \
      TMP_DIR=${TMPDIR}
    """
}

