#! /usr/bin/env nextflow
/*************************************************
 * This module wraps any scripts in bin in processes.
 * Eventually the scripts will be split into smaller processes for parallization
 *************************************************/

/*
 * gatk0_index.sh help statement
 */
 
process gatk0_index_help {

  output: path 'gatk0_index_help.txt'

  script:
  """
  bash gatk0_index.sh > gatk0_index_help.txt
  """
}

/*
 * bash gatk0_index.sh <genome.fasta> <name_str>
 */
 
process gatk0_index {
// Note: labels will link to configs for loading modules, containers, other run time info
// this is a placeholder
//   label "samtools"
//   label "picard"
//   label "bwa"
//   label "bedtools2"
//   label "bioawk"

  label "gatk0_index"

  input:
  path genome_fasta
  
  output:
  path 'gatk0_index.out'

  script:
  """
  echo "I can see $genome_fasta" > gatk0_index.out
  # bash gatk0_index.sh "$genome_fasta" "genome"
  """
}
