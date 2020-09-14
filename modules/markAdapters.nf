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
    tag "$read_bam"
    label 'picard'
    publishDir "${params.outdir}/markAdapters", mode: 'copy'

    input:
    path read_bam
    //val TMPDIR

    output:
    path "${read_bam.baseName}*.bam", emit: read_marked
    path "${read_bam.baseName}*.txt"
//    path "readname_markilluminaadapters_metrics.txt"

    script:
    """
    picard MarkIlluminaAdapters \
      I=$read_bam \
      O=${read_bam.baseName}_markilluminaadapters.bam \
      M=${read_bam.baseName}_markilluminaadapters_metrics.txt
    """
}

/* Scrap
TMP_DIR=$TMPDIR
*/
