#! /usr/bin/env nextflow

process gatk {
  label 'gatk'
  publishDir 'gatk'

  input:
  path genome_fasta
  path reads_bams
  path coord_list
  
  output:
  path "*.vcf", emit: vcf
  path "*.idx", emit: idx

  script:
  """
  #! /usr/bin/env bash
  gatk --java-options "-Xmx80g -XX:+UseParallelGC" HaplotypeCaller \
  -R ${genome_fasta} \
  -I reads_bams \
  -L chr1:1-999999 \
  --output ${coord_list.simpleName}.vcf;
  """
}