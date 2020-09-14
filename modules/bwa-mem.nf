#! /usr/bin/env nextflow

bwa_container = 'kathrinklee/bwa'

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
    tag "$genome_fasta, $readname"
    label 'bwa'
    publishDir "${params.outdir}/bwa-mem", mode: 'copy'

    input:
    path genome_fasta
    path genome_index
    each readname_fq

    output:
    path "$genome_fasta"
    path "${readname_fq.baseName}*.bam"

    script:
    """
    bwa mem \
      -M \
      -t 15 \
      -p $genome_fasta \
    ${readname_fq} |\
   samtools view -buS - > ${readname_fq.baseName}_bwa_mem.bam
    """
}
