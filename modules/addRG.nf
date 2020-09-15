#! /usr/bin/env nextflow

picard_container = ''

process AddOrReplaceReadGroups_help {
  label 'picard'

  container = "$picard_container"

  output: path 'picard-AddOrReplaceReadGroups.txt'

  """
  picard AddOrReplaceReadGroups > picard-AddOrReplaceReadGroups.txt
  """
}

process AddOrReplaceReadGroups_run {
    tag "$readname"
    label 'picard'
    publishDir "${params.outdir}/AddOrReplaceReadGroups", mode: 'copy'

    input:
    tuple val(readname), path(readpairs)
    val(RGID)
    val(RGLB)
    val(RGPL)
    val(RGPU)
    val(RGSM)
    val(TMPDIR)

    output:
    path "${readname}_final.bam"

    script:
    """
    #! /usr/bin/env bash
    picard AddOrReplaceReadGroups \
    INPUT=${readname}_MarkDuplicates.bam \
	OUTPUT=${readname}_final.bam \
    RGID=${RGID} \
    RGLB=${RGLB} \
    RGPL=${RGPL} \
    RGPU=${RGPU} \
    RGSM=${RGSM} \
    CREATE_INDEX=true \
    TMP_DIR=$TMPDIR
    """
}

/* Scrap here


*/
