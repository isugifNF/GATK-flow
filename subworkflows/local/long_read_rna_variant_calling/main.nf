#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { longcallR } from '../../../modules/local/longcallR.nf'
include { pbmm2_index; pbmm2_align; gatk_HaplotypeCaller as gatk_HaplotypeCaller_LR_RNA } from '../../../modules/local/LongReadseq.nf'
include { CreateSequenceDictionary; 
          samtools_faidx; 
          MarkDuplicates;
          bedtools_makewindows; } from '../../../modules/local/GATK.nf'
include { SplitNCigarReads } from '../../../modules/local/RNAseq.nf'

workflow LONGREAD_RNA_VARIANT_CALLING {
  take:
  genome_ch
  reads_ch

  main:
  // == Since one sample may be run on multiple lanes
  i = 1

  // == Prepare mapped read files
  cleanreads_ch = reads_ch
    | map { n -> ["${i++}_"+n.baseName, n] }

  genome_ch 
    | (CreateSequenceDictionary & samtools_faidx )

  // == Alignment with pbmm2
  genome_ch
    | pbmm2_index
    | combine(cleanreads_ch)
    | pbmm2_align

  // == GATK Best Practices Preprocessing
  // MarkDuplicates expects tuple val(i_readname), path(merge_bam), path(merge_bai)
  processed_ch = pbmm2_align.out
    | MarkDuplicates
    | combine(genome_ch)
    | combine(samtools_faidx.out)
    | combine(CreateSequenceDictionary.out)
    | SplitNCigarReads

  // == Variant Calling
  if (params.lr_rna_caller == 'longcallR') {
    // Run longcallR
    vcf_out_ch = processed_ch
      | combine(genome_ch)
      | longcallR
  } else {
    // Run GATK HaplotypeCaller
    // Split windows for GATK
    windows_ch = samtools_faidx.out
      | bedtools_makewindows
      | splitText(){it.trim()}

    vcf_out_ch = processed_ch
      | combine(windows_ch)
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | gatk_HaplotypeCaller_LR_RNA
  }

  emit:
  vcf = vcf_out_ch
}
