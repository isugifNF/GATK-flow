#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/* import modules */
include { gatk0_index; gatk1_preprocess } from './modules/wrap_bin.nf'
include { fastqc } from './modules/fastqc.nf'

include { sortSeq_run }      from './modules/sortSeq.nf'
include { createSeqDict_run} from './modules/createSeqDict.nf'
include { bwa_index }        from './modules/bwa.nf'
include { seqLength_run }    from './modules/seqLength.nf'
include { faidx_run }        from './modules/faidx.nf'
include { bedtools_coords }  from './modules/makeIntervals.nf'
include { fastqToSAM_run }   from './modules/fastqToSAM.nf'

/* define workflow */

workflow {
  //==== (Step gatk-1) quality check reads
  //channel.fromPath(params.reads, checkIfExists:true) | fastqc

  //==== (Step gatk0) index genome for faster alignment
  // channel.fromPath(params.genome) | gatk0_index
  channel.fromPath(params.genome, checkIfExists:true) | sortSeq_run |\
     (createSeqDict_run & bwa_index & seqLength_run & faidx_run)

  seqLength_run.out | bedtools_coords

  //==== (Step gatk1) prepare
  channel.fromFilePairs(params.reads, checkIfExists:true) | fastqToSAM_run

//channel.fromFilePairs(params.reads, checkIfExists:true) |
//    flatten | buffer(size:2, skip:1) |
//    flatten | fastqc

  //.view { it.key + ': ' + it.value }
  //channel.fromFilePairs(params.reads) | fromFilePairs_toFiles | flatten | fastqc
  //paired_files | view

}

// =========== Scrap
// gatk1 preprocess
//readpairs_ch = channel.fromFilePairs(params.reads)
//gatk1_preprocess(
//    gatk0_index.out.indexed_genome,
//    channel.fromFilePairs(params.reads)
//)
//gatk1_preprocess(
//    gatk0_index.out.indexed_genome,
//  channel.fromPath(params.reads).take(2).flatten()
//                 )
//gatk2_preprocess_help()
//gatk3_cmdsgen_help()
//gatk4_filter_help()

// Debug check: make sure files are passed in
//channel.fromPath(params.genome) | view
//channel.fromPath(params.reads).take(2) | view
//  channel.fromFilePairs(params.reads) | println
