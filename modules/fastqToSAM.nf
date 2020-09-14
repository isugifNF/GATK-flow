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
    tag "$readname"
    label 'picard'
    publishDir "${params.outdir}/fastqToSAM", mode: 'copy'

    input:
    tuple val(readname), path(readpairs)
//    val read1
//    val read2
//    val readgroup
//    val readname
//    val platform
//    val rundate
//    var center

    output:
    path "${readname}_fastqtosam.bam"

    script:
    """
    #! /usr/bin/env bash
    picard FastqToSam \
      FASTQ=${readpairs.get(0)} \
      FASTQ2=${readpairs.get(1)} \
      OUTPUT=${readname}_fastqtosam.bam \
      SAMPLE_NAME=${readname}
    """
}

/* Scrap here
#platform_unit=`gunzip -c \$read1 | head -n 1 | cut -f 3 -d ":"`
#center="ISU"
#rundate=\$(date '+%Y-%m-%d %H:%M:%S' |sed 's/ /T/g')
  #    READ_GROUP_NAME={readname} \
  #    LIBRARY_NAME={readname}_lib \
  #   PLATFORM_UNIT={platform_unit} \
  #    PLATFORM=ILLUMINA \
  #    SEQUENCING_CENTER={center} \
  #    RUN_DATE={rundate}
*/
