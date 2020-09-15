#! /usr/bin/env nextflow

picard_container = ''

process MarkDuplicates_help {
  label 'picard'

  container = "$picard_container"

  output: path 'picard-MarkDuplicates.txt'

  """
  picard MarkDuplicates > picard-MarkDuplicates.txt
  """
}

process MarkDuplicates_run {
    tag "$readname"
    label 'picard'
    publishDir "${params.outdir}/MarkDuplicates", mode: 'copy'

    input:
    tuple val(readname), path(readpairs)

    output:
    path "${readname}_MarkDuplicates.bam"

    script:
    """
    #! /usr/bin/env bash
    picard MarkDuplicates \
    INPUT=${readname}_MergeBamAlignment.bam \
    OUTPUT=${readname}_MarkDuplicates.bam \
    METRICS_FILE=${readname}_markduplicates_metrics.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    CREATE_INDEX=true
    """
}

/* Scrap here



*/
