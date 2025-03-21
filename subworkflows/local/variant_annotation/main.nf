#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { snpEff_Build; snpEff_Annotate; } from '../../../modules/local/snpEff.nf'

workflow SNP_EFF_ANNOTATION {
  take:
  genome_ch
  gff_ch
  vcf_ch

  main:

//  annotated_vcf_ch = genome_ch
//  | combine(gff_ch)
//  | snpEff_Build    // Building snpEff database
//  | combine(genome_ch)
//  | combine(vcf_ch)
//  | snpEff_Annotate // Annotating variants using snpEff

  snpEff_db_ch = genome_ch 
    | combine(gff_ch)
    | snpEff_Build

  // Create combinations for annotation 
  annotation_input_ch = snpEff_db_ch
    | map { db -> tuple(db, db.name + ".fa") }
    | combine(vcf_ch)
    | map { db, pattern, vcf -> tuple(db, file(pattern), vcf) }

  // Annotate variants
  annotated_vcf_ch = annotation_input_ch
    | snpEff_Annotate

  emit:
  annotated_vcf = annotated_vcf_ch
}
