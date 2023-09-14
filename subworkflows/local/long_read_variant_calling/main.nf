#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { CreateSequenceDictionary;
          samtools_faidx;
          bedtools_makewindows;
          CombineGVCFs;
          GenotypeGVCFs;
          SortVcf;
          calc_DPvalue;
          VariantFiltration;
          keep_only_pass; } from '../../../modules/local/GATK.nf'

include { pbmm2_index; 
          pbmm2_align;
          gatk_HaplotypeCaller as gatk_HaplotypeCaller_LongRead;
          CombineGVCFs as CombineGVCFs_LongRead; } from '../../../modules/local/longReadseq.nf'

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
    | map { n -> ["${i++}_"+n.basename, n] }

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
      | combine(calc_DPvalue.out.map{n-> n.replaceAll("\n","")})
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | VariantFiltration
      | keep_only_pass

  emit:
  vcf = vcf_out_ch
}
