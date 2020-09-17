#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/* import modules */
// include { get_test_data } from './modules/wrap_bin.nf'
include { fastqc } from './modules/fastqc.nf'

include {
  sortSeq_run
  seqLength_run
} from './modules/bioawk.nf'

include {
  createSeqDict_run
  fastqToSAM_run
  markAdapters_run
  SamToFastq_run
} from './modules/picard.nf'
include {
  bwa_index
  bwa_mem_run
} from './modules/bwa.nf'

include { faidx_run }        from './modules/faidx.nf'
include { bedtools_coords }  from './modules/makeIntervals.nf'
include { MergeBamAlignment_run } from './modules/mergeAlignment.nf'

/* define workflow */

workflow {
  //==== (Run once to fetch test data)
  //  get_test_data()

  //==== (Step gatk-1) quality check reads
  //channel.fromPath(params.reads, checkIfExists:true) | fastqc

  //==== (Step gatk0) index genome for faster alignment
  // channel.fromPath(params.genome) | gatk0_index
  channel.fromPath(params.genome, checkIfExists:true) |\
     sortSeq_run |\
     (createSeqDict_run & bwa_index & seqLength_run & faidx_run)

  seqLength_run.out | bedtools_coords

  //==== (Step gatk1) prepare
  channel.fromFilePairs(params.reads, checkIfExists:true).take(3) |\
    fastqToSAM_run |\
    markAdapters_run

  markAdapters_run.out.read_marked | SamToFastq_run


//  createSeqDict_run.out.view()

  bwa_mem_run(
    bwa_index.out,
    SamToFastq_run.out
  )

  //bwa_index.out.mix().view()
  //createSeqDict_run.out | view

  genome_ch = bwa_index.out.mix().collate(2,2).combine(createSeqDict_run.out)

  /* This makes sure mapped and unmapped read names match, or it will combine randomly */
  read_ch1 = fastqToSAM_run.out
    .flatMap {n -> [n.baseName, n] }
    .collate(2,2)

  read_ch2 = bwa_mem_run.out
    .flatMap { n -> [ n.baseName.replaceFirst("_markilluminaadapters_interleaved_bwa_mem", ""), n] }
    .collate(2,2)

  read_merge_ch3 = read_ch1.join(read_ch2)

  read_merge_ch3.combine(genome_ch) | view
  /*
  read_merge_ch3
    .flatMap{ n -> [n, bwa_index.out.value()]} | view
    */

/*
  MergeBamAlignment_run (
    bwa_index.out,
    createSeqDict_run.out,
    read_merge_ch3
  )
*/

}
