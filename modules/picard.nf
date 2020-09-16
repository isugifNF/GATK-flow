#! /usr/bin/env nextflow

//picard_container = ''

process createSeqDict_help {
  label 'picard'

  output: path 'picard-seqCreateDict.txt'

  """
  picard CreateSequenceDictionary > picard-seqCreateDict.txt
  """
}

process createSeqDict_run {
    tag "$sorted_ref"
    label 'picard'
    publishDir "${params.outdir}/createSeqDict", mode: 'copy'

    input:
    path sorted_ref

    output:
    path "${sorted_ref.baseName}.dict"

    script:
    """
    picard CreateSequenceDictionary \
      REFERENCE=${sorted_ref} \
      OUTPUT=${sorted_ref.baseName}.dict
    """
}

process fastqToSAM_help {
  label 'picard'

//  container = "$picard_container"

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
#READ_GROUP_NAME={readname} \
#LIBRARY_NAME={readname}_lib \
#PLATFORM_UNIT={platform_unit} \
#PLATFORM=ILLUMINA \
#SEQUENCING_CENTER={center} \
#RUN_DATE={rundate}
*/

process markAdapters_help {
  label 'picard'

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

process SamToFastq_help {
  label 'picard'

  container = "$picard_container"

  output: path 'picard-SamToFastq.txt'

  """
  picard SamToFastq > picard-SamToFastq.txt
  """
}

process SamToFastq_run {
    tag "$read_marked"
    label 'picard'
    publishDir "${params.outdir}/SamToFastq", mode: 'copy'

    input:
    path read_marked

    output:
    path "${read_marked.baseName}_samtofastq_interleaved.fq"

    script:
    """
    picard SamToFastq \
      I=${read_marked} \
      FASTQ=${read_marked.baseName}_samtofastq_interleaved.fq \
      CLIPPING_ATTRIBUTE=XT \
      CLIPPING_ACTION=2 \
      INTERLEAVE=true \
      NON_PF=true
    """
}
