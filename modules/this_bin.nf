#! /usr/bin/env nextflow


process say_hi {

  output: path 'say_hi_out.txt'

  script:
  """
  bash say_hi.sh > say_hi_out.txt
  """
}


/* this is a placeholder process, we can split this script later */

process gatk0_index {
// Note: labels will link to configs for loading modules, containers, other run time info
// this is a placeholder
//   label "samtools"
//   label "picard"
//   label "bwa"
//   label "bedtools2"
//   label "bioawk"

//  input:
  output: path 'gatk0_index_out.txt'

  script:
  """
  gatk0_index.sh > gatk0_index_out.txt
  """

}
