#! /usr/bin/env nextflow

process AddOrReplaceReadGroups_help {
  label 'picard'

  output:
  path 'picard-AddOrReplaceReadGroups.txt'

  """
  picard AddOrReplaceReadGroups > picard-AddOrReplaceReadGroups.txt
  """
}

process AddOrReplaceReadGroups_run {
    tag "$marked_bam.fileName"
    label 'picard'
    publishDir "${params.outdir}/AddOrReplaceReadGroups", mode: 'copy'

    input:
    path marked_bam

    output:
    path "${marked_bam.simpleName}_final.bam"

    script:
    """
    #! /usr/bin/env bash
    picard AddOrReplaceReadGroups \
    INPUT=${marked_bam} \
    OUTPUT=${marked_bam.simpleName}_final.bam \
    RGID=4 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=20 \
    CREATE_INDEX=true
    """
}

/* Scrap here
/*
    tuple val(readname), path(readpairs)
    val(RGID)
    val(RGLB)
    val(RGPL)
    val(RGPU)
    val(RGSM)
    val(TMPDIR)

java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
      O=output.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
*/
