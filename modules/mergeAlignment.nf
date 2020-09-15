#! /usr/bin/env nextflow

picard_container = ''

process MergeBamAlignment_help {
  label 'picard'

  container = "$picard_container"

  output: path 'picard-MergeBamAlignment.txt'

  """
  picard MergeBamAlignment > picard-MergeBamAlignment.txt
  """
}

process MergeBamAlignment_run {
    tag "$readname"
    label 'picard'
    publishDir "${params.outdir}/MergeBamAlignment", mode: 'copy'

    input:
    tuple val(readname), path(readpairs)
    val genome

    output:
    path "${readname}_MergeBamAlignment.bam"

    script:
    """
    #! /usr/bin/env bash
    picard MergeBamAlignment \
    R=$genome \
    UNMAPPED_BAM=${readname}_fastqtosam.bam \
    ALIGNED_BAM=${readname}_bwa_mem.bam \
    O=${readname}_MergeBamAlignment.bam \
    CREATE_INDEX=true \
    ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \
    CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    ATTRIBUTES_TO_RETAIN=XS
    """
}

/* Scrap here


*/
