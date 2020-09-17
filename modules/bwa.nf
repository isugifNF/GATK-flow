#! /usr/bin/env nextflow

bwa_container = 'kathrinklee/bwa'

process bwa_help {
  label 'bwa'

  container = "$bwa_container"

  output: path 'bwa_out.txt'

  """
  echo "bwa" > bwa_out.txt
  bwa &>> bwa_out.txt
  """
}

process bwa_index {
    tag "$sorted_ref"
//    label 'process_medium'
    label 'bwa'
    publishDir "${params.outdir}/bwa", mode: 'copy'

//    when:
//    !params.skipQC && !params.skipFastQC

    input:
//    tuple val(name), file(reads) //from raw_reads_fastqc
    path(sorted_ref)

    output:
    path "$sorted_ref"
    path "$sorted_ref*"

    script:
    """
    bwa index $sorted_ref
    """
}

process bwa_mem_help {
  label 'bwa'

  container = "$bwa_container"

  output: path 'bwa-mem_out.txt'

  """
  echo "bwa mem" > bwa-mem_out.txt
  bwa mem &>> bwa-mem_out.txt
  """
}

process bwa_mem_run {
    tag "$genome_fasta, $genome_index, $readname_fq.simpleName"
    label 'bwa'
    label 'samtools'
    publishDir "${params.outdir}/bwa_mem", mode: 'copy'

    input:
    path genome_fasta
    path genome_index
    each readname_fq

    output:
    path "$genome_fasta"
    path "${readname_fq.simpleName}*.bam"

    script:
    """
    bwa mem \
      -M \
      -t 15 \
      -p $genome_fasta \
    ${readname_fq} |\
   samtools view -buS - > ${readname_fq.simpleName}_bwa_mem.bam
    """
}
