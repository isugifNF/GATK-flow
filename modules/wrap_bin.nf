#! /usr/bin/env nextflow
/*************************************************
 * This module wraps any scripts in bin in processes.
 * Eventually the scripts will be split into smaller processes for parallization
 *************************************************/

/* Fetches test-data folder, really this is a note to myself"
// process get_test_data {
//   publishDir "$params.outdir/test-data", mode:'copy'
// 
//   output:
//   path "test-data"
//   path "test-data/*"
//   
//   script:
//   """
//   # /usr/bin/env bash
//   wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
//   tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
//   """
// }


/*
 * gatk0_index.sh help statement
 */

process gatk0_index_help {
  publishDir "$params.outdir/wrap_bin"

  output: path 'gatk0_index_help.txt'

  script:
  """
  #! /usr/bin/env bash
  bash gatk0_index.sh > gatk0_index_help.txt
  """
}

/*
 * bash gatk0_index.sh <genome.fasta> <name_str>
 */

process gatk0_index {
  tag "$genome_fasta"

//  label 'samtools'
//  label 'picard'
//  label 'bwa'
//  label 'bedtools2'
//  label 'bioawk'

  label "gatk0_index"
  publishDir "${params.outdir}/gatk0_index", mode: 'copy'

  input:
  path genome_fasta

  output:
  path "$genome_fasta", emit: indexed_genome
  path 'gatk0_index.out'
  path "${genome_fasta.baseName}*"


  script:
  """
  echo "I can see $genome_fasta" > gatk0_index.out
  bash gatk0_index.sh "$genome_fasta" "${genome_fasta.baseName}"
  """
}

process gatk1_preprocess {
  tag "$genome_fasta, $name, $reads"

//  label 'samtools'
//  label 'picard'
//  label 'bwa'

  label "gatk1_index"

  publishDir "${params.outdir}/gatk1_preprocess", mode: 'copy'

  input:
  path genome_fasta
  tuple val(name), path(reads)

  output:
  path "$name*"

  script:
  """
  bash gatk1_preprocess.sh $genome_fasta $name $reads
  """
}

/*
 * gatk3_cmdsgen.sh
 */

process gatk3_cmdsgen_help {
  output: path 'gatk3_cmdsgen_help.txt'

  script:
  """
  bash gatk3_cmdsgen.sh > gatk3_cmdsgen_help.txt
  """
}

/*
 * gatk4_filter.sh
 */

process gatk4_filter_help {
  output: path 'gatk4_filter_help.txt'

  script:
  """
  bash gatk4_filter.sh > gatk4_filter_help.txt
  """
}
