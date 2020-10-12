

process merge_vcf {
  publishDir "${params.outdir}/vcftools"

  input:
  tuple path(vcfs)

  output:
  path "first-round_merged.vcf"

  script:
  """
  #! /usr/bin/env bash
  cat ${vcfs.get(0)} | grep "^#" > first-round_merged.vcf
  cat $vcfs | grep -v "^#" >> first-round_merged.vcf
  """
}

process vcftools_snp_only {
  tag "${merged_vcf.fileName}"
  label 'vcftools'
  publishDir "${params.outdir}/vcftools"

  input:
  path merged_vcf

  output:
  path "${merged_vcf.simpleName}_snps-only.*"

  script:
  """
  #! /usr/bin/env bash
  $vcftools_app \
    --vcf $merged_vcf \
      --remove-indels \
      --recode \
      --recode-INFO-all \
      --out ${merged_vcf.simpleName}_snps-only
  """
}

process run_SortVCF {
  label 'picard'
  publishDir "$params.outdir/picard"

  input:
  tuple path(vcf), path(dict)

  output:
  path ("$vcf.simpleName.sorted*")

  script:
  """
  #! /usr/bin/env bash
  $picard_app SortVcf \
    INPUT=$vcf \
    SEQUENCE_DICTIONARY=$dict \
    CREATE_INDEX=true \
    OUTPUT=${vcf.simpleName}.sorted.vcf
  """
}

// java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar

process calc_DPvalue {
  tag "$sorted_vcf.fileName"
  label 'datamash'

  input:
  path(sorted_vcf)

  output:
  stdout()

  script:
  """
  #! /usr/bin/env bash
  grep -v "^#" $sorted_vcf | cut -f 8 | grep -oe ";DP=.*" | cut -f 2 -d ';' | cut -f 2 -d "=" > dp.txt
  cat dp.txt | datamash mean 1 sstdev 1 > dp.stats
  cat dp.stats | awk '{print \$1+5*\$2}'
  """
}

process gatk_VariantFiltration {
  input:
  tuple path(sorted_snp_vcf), val(dp), path(genome_fasta), path(genome_index), path(genome_fai), path(genome_dict)

  output:
  path("${sorted_snp_vcf.simpleName}.marked.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app VariantFiltration \
    --reference $genome_fasta \
    --variant $sorted_snp_vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > $dp" \
    --filterName "FAIL" \
    --output ${sorted_snp_vcf.simpleName}.marked.vcf
  """
}

// java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar

process keep_only_pass {
  input:
  path(snp_marked_vcf)

  output:
  path("${snp_marked_vcf.simpleName}_snp-only.pass_only_vcf")

  script:
  """
  #! /usr/bin/env bash
  grep -v "^#" $snp_marked_vcf | awk '\$7=="PASS"' > body.pass
  grep "^#" $snp_marked_vcf > head.pass
  cat head.pass body.pass > ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  """
}
