#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/* import modules */
include { gatk0_index_help; gatk0_index; gatk2_preprocess_help; gatk3_cmdsgen_help; gatk4_filter_help } from './modules/wrap_bin.nf'
include { fastqc } from './modules/fastqc.nf'

/* define workflow */

workflow {
  // check quality of fastq files
  channel.fromPath(params.reads) | fastqc

  // index genome
  gatk0_index_help()
  channel.fromPath(params.genome) | gatk0_index
  gatk2_preprocess_help()
  gatk3_cmdsgen_help()
  gatk4_filter_help()

  // Debug check: make sure files are passed in
  channel.fromPath(params.genome) | view
  channel.fromPath(params.reads).take(3) | view
//  channel.fromFilePairs(params.reads) | println
}

