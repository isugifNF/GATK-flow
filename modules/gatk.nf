#! /usr/bin/env nextflow

/* todo: parallize the while loop */

process gatk {
  label 'gatk'
  publishDir 'gatk'

  input:
  tuple path genome_fasta, path genome_indexes
  tuple path reads_bams, path reads_bai
  path coord_list
  
  output:
  path "*.vcf", emit: vcf
  path "*.idx", emit: idx

  script:
  """
  #! /usr/bin/env bash
  
  while read LINE; do
    gatk --java-options "-Xmx80g -XX:+UseParallelGC" HaplotypeCaller \
    -R ${genome_fasta} \
    -I $reads_bams \
    -L \${LINE} \
    --output ${reads_bam.simpleName}_\${LINE}.vcf
  done < $coord_list
  """
}
