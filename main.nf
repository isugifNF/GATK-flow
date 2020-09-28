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
include { MergeBamAlignment_run; test_run } from './modules/mergeAlignment.nf'
include { MarkDuplicates_run } from './modules/mergeDuplicates.nf'
include { AddOrReplaceReadGroups_run } from './modules/addRG.nf'
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

  // Single genome channel

  genome_ch = bwa_index.out.genome_fasta
    .combine(bwa_index.out.genome_index.toList())
    .combine(createSeqDict_run.out)
    .combine(faidx_run.out)

  // Multiple read channel
  /* This makes sure mapped and unmapped read names match, or it will combine randomly */

  read_unmapped = fastqToSAM_run.out
    .flatMap {
      // Create a list [BiosampleX_string, unmapped_bam_file]
      n -> [n.simpleName, n]
    }
    .collate(2)

  read_mapped = bwa_mem_run.out
    .flatMap {
      // Create a list [BiosampleX_string, mapped_file]
      n -> [ n.simpleName.replaceFirst("_markilluminaadapters_interleaved_bwa_mem", ""), n]
    }
    .collate(2)


  reads_ch = read_unmapped
    .join(read_mapped)
    .combine(genome_ch)

  reads_ch | MergeBamAlignment_run 
  MergeBamAlignment_run.out.bam | MarkDuplicates_run
  MarkDuplicates_run.out.bam | AddOrReplaceReadGroups_run

  AddOrReplaceReadGroups_run.out.flatten() | view

}
