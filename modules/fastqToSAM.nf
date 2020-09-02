#! /usr/bin/env nextflow

picard_container = ''

process fastqToSAM_help {
  label 'picard'

  container = "$picard_container"
  
  output: path 'picard-fastqToSAM.txt'

  """
  picard FastqToSam > picard-fastqToSAM.txt
  """
}

process fastqToSAM_run {
    tag "$name"
    label 'picard'
    publishDir "${params.outdir}/fastqToSAM", mode: 'copy'

    input:
    val read1
    val read2
    val readgroup
    val readname
    val platform
    val rundate
    var center 

    output:
    path "readname_fastqtosam.bam"

    script:
    """
    picard FastqToSam \
      FASTQ=${read1} \
      FASTQ2=${read2} \
      OUTPUT=${readname}_fastqtosam.bam \
      READ_GROUP_NAME=${readgroup} \
      SAMPLE_NAME=${readname} \
      LIBRARY_NAME=${readname}_lib \
      PLATFORM_UNIT=${platform} \
      PLATFORM=illumina \
      SEQUENCING_CENTER=${center} \
      RUN_DATE=${rundate}
    """
}
