#! /usr/bin/env nextflow

bwa_container = 'kathrinklee/bwa'

process bwa-mem_help {
  label 'bwa'

  container = "$bwa_container"
  
  output: path 'bwa-mem_out.txt'

  """
  echo "bwa mem" > bwa-mem_out.txt
  bwa mem &>> bwa-mem_out.txt
  """
}

process bwa-mem_run {
    tag "$name"
    label 'bwa'
    publishDir "${params.outdir}/bwa-mem", mode: 'copy'

    input:
    path(genome_fasta)
    val readname

    output:
    path "$genome_fasta"
    path "*.bam"

    script:
    """
    bwa mem \
      -M \
      -t 15 \
      -p $genome_fasta \
    ${readname}_samtofastq_interleaved.fq |\
   samtools view -buS - > ${readname}_bwa_mem.bam 
    """
}
