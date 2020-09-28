#! /usr/bin/env nextflow

picard_container = ''

process MarkDuplicates_help {
  label 'picard'

  output: path 'picard-MarkDuplicates.txt'

  """
  picard MarkDuplicates > picard-MarkDuplicates.txt
  """
}

process MarkDuplicates_run {
    tag "${merged_bam.fileName}"
    label 'picard'
    publishDir "${params.outdir}/MarkDuplicates", mode: 'copy'

    input:
    path merged_bam

    output:
    path "${merged_bam.simpleName}*"

    script:
    """
    #! /usr/bin/env bash
    picard MarkDuplicates \
    INPUT=${merged_bam} \
    OUTPUT=${merged_bam.simpleName}_MarkDuplicates.bam \
    METRICS_FILE=${merged_bam.simpleName}_markduplicates_metrics.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    CREATE_INDEX=true
    """
}

/* Scrap here



*/
