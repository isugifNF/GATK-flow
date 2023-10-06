#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { CreateSequenceDictionary;
          samtools_faidx;
          calc_DPvalue;
          keep_only_pass; } from '../../../modules/local/GATK.nf'

include { pbmm2_index; 
          pbmm2_align;
          gatk_HaplotypeCaller as gatk_HaplotypeCaller_LongRead;
          CombineGVCFs as CombineGVCFs_LongRead;
          GenotypeGVCFs as GenotypeGVCFs_LongRead;
          calc_DPvalue as calc_DPvalue_LongRead;
          VariantFiltration as VariantFiltration_LongRead } from '../../../modules/local/LongReadseq.nf'

include { MarkDuplicates as MarkDuplicates_RNA; } from '../../../modules/local/RNAseq.nf'



workflow LONGREAD_VARIANT_CALLING {
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

  mapped_ch = genome_ch
    | pbmm2_index
    | combine(cleanreads_ch)
    | pbmm2_align
    | MarkDuplicates_RNA
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)
    | gatk_HaplotypeCaller_LongRead 
    | collect
    | map { n -> [n] }
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)
    | CombineGVCFs_LongRead
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)
    | GenotypeGVCFs_LongRead

  CombineGVCFs_LongRead.out
    | calc_DPvalue_LongRead

  vcf_out_ch = GenotypeGVCFs_LongRead.out
    | combine(calc_DPvalue_LongRead.out.map{n-> n.replaceAll("\n","")})
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)
    | VariantFiltration_LongRead
    | keep_only_pass

  emit:
  vcf = vcf_out_ch
}
