#! /usr/bin/env nextflow
nextflow.enable.dsl=2

process longcallR {
  tag "$readname"
  label 'longcallR'
  container params.longcallR_container
  publishDir "${params.outdir}/04_longcallR", mode: 'copy'

  input:
  tuple val(readname), path(bam), path(bai), path(genome_fasta)

  output:
  tuple val(readname), path("${readname}_longcallR.vcf")

  script:
  """
  # longcallR typically expects a BAM file and a reference
  # Adjust command based on specific longcallR documentation
  # Assuming standard usage: longcallR -b input.bam -f ref.fa -o output_prefix
  
  longcallR \
    -b ${bam} \
    -f ${genome_fasta} \
    -o ${readname}_longcallR
    
  # If the output is not exactly .vcf, rename or adjust capture
  if [ -f "${readname}_longcallR.vcf.gz" ]; then
      gunzip "${readname}_longcallR.vcf.gz"
  fi
  """

  stub:
  """
  touch ${readname}_longcallR.vcf
  """
}
