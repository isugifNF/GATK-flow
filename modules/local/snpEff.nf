#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process snpEff_Build {
  tag "${genome_fasta.simpleName}"
  label 'snpeff'
  publishDir "${params.outdir}/snpEff_Build"
  
  input:
  tuple path(genome_fasta), path(gff_file)

  output:
  path("data/${genome_fasta.simpleName}")

  script:
  """
  mkdir -p data/${genome_fasta.simpleName}
  cp $gff_file data/${genome_fasta.simpleName}/genes.gff
  cp $genome_fasta data/${genome_fasta.simpleName}/sequences.fa
  java -jar /opt/snpEff/snpEff.jar build -configOption ${genome_fasta.simpleName}.genome=${genome_fasta.simpleName} -dataDir `pwd`/data/ -gff3 -v -noCheckCds -noCheckProtein ${genome_fasta.simpleName}
  """
}

process snpEff_Annotate {
  tag "${vcf_file.simpleName}"
  label 'snpeff'
  publishDir "${params.outdir}/snpEff_Annotate"

  input:
  tuple path(snpeff_db), path(genome_fasta), path(vcf_file)

  output:
  tuple path("${vcf_file.simpleName}.eff.vcf"), path("${vcf_file.simpleName}_summary.html"), path("${vcf_file.simpleName}_genes.txt")

  script:
  """
  java -jar /opt/snpEff/snpEff.jar eff \
    -configOption ${snpeff_db}.genome=${snpeff_db} \
    -v \
    -dataDir ${snpeff_db} \
    $vcf_file \
    > ${vcf_file.simpleName}.eff.vcf

  mv snpEff_summary.html ${vcf_file.simpleName}_summary.html
  mv snpEff_genes.txt ${vcf_file.simpleName}_genes.txt
  """
}
