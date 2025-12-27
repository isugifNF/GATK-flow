#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { longcallR } from '../../../modules/local/longcallR.nf'
include { pbmm2_index; pbmm2_align } from '../../../modules/local/LongReadseq.nf'
include { CreateSequenceDictionary; samtools_faidx } from '../../../modules/local/GATK.nf'

workflow LONGREAD_RNA_VARIANT_CALLING {
  take:
  genome_ch
  reads_ch

  main:
  // == Since one sample may be run on multiple lanes
  i = 1

  // == Prepare mapped and unmapped read files
  cleanreads_ch = reads_ch
    | map { n -> ["${i++}_"+n.baseName, n] }

  genome_ch 
    | (CreateSequenceDictionary & samtools_faidx )

  // Reuse pbmm2 for alignment (Note: Ensure pbmm2 parameters are suitable for RNA/Iso-Seq if needed)
  // Alternatively, implement minimap2 if pbmm2 is insufficient for splice alignment
  vcf_ch = genome_ch
    | pbmm2_index
    | combine(cleanreads_ch)
    | pbmm2_align
    | map { n -> [n[0], n[1], n[2], n[3]] } // Adjust tuple if needed, pbmm2_align output: [val(readname), path(bam), path(bai)]
    // pbmm2_align output is tuple val("$readname"), path("${readname}_mapped.bam"), path("${readname}_mapped.bam.bai")
    // We need to combine this with genome_fasta for longcallR
    
  // Re-map to get [readname, bam, bai]
  aligned_ch = pbmm2_align.out
  
  // Combine with genome
  longcallR_input_ch = aligned_ch
    | combine(genome_ch)
    
  // Run longcallR
  longcallR_input_ch
    | longcallR

  emit:
  vcf = longcallR.out
}
