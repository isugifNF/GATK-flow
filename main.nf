#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/* import modules */
include { gatk0_index_help; gatk0_index } from './modules/wrap_bin.nf'

/* define workflow */

workflow {
  gatk0_index_help()
  channel.fromPath(params.genome) | gatk0_index

  /* Make sure files are passed in */
  channel.fromPath(params.genome) | view
  channel.fromPath(params.reads).take(3) | view
//  channel.fromFilePairs(params.reads) | println
}

