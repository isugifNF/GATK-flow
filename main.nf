#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/* import modules */
include { say_hi; gatk0_index } from './modules/this_bin.nf'

/* define workflow */

workflow {
  say_hi()
//  gatk0_index()
}

