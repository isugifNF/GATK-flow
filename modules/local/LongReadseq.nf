#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process pbmm2_index {
  tag "${genome_fasta.simpleName}"
  label 'pbmm2'
  publishDir "${params.outdir}/02_MapReads"

  input:
  path(genome_fasta)

  output: // [genome.fasta, [genome_index files]]
  tuple path("$genome_fasta"), path("${genome_fasta}*")

  script:
  """
  #! /usr/bin/env bash
  $pbmm2_app index $genome_fasta ${genome_fasta}.mmi
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta}.mmi
  """
}

process pbmm2_align {
  tag "$readname"
  label 'pbmm2'
  publishDir "${params.outdir}/02_MapReads"

  input:
  tuple path(genome_fasta), path(genome_index), val(readname), path(reads)

  output: // reads_mapped_2_genome.bam
  tuple val("$readname"), path("${readname}_mapped.bam"), path("${readname}_mapped.bam.bai")

  script:
  """
  #! /usr/bin/env bash
  $pbmm2_app align ${genome_index} \
    ${reads} \
    ${readname}_mapped.bam \
    --preset CCS \
    --sort \
    --rg '@RG\tID:ZP\tSM:zapped'
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${readname}_mapped.bam
  touch ${readname}_mapped.bam.bai
  """
}

process gatk_HaplotypeCaller {
  tag "$bam"
  label 'gatk_long'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'
  errorStrategy 'retry'
  maxRetries 3

  input:  // [window, reads files ..., genome files ...]
  tuple val(read_name), path(bam), path(bai), val(window), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  tuple path("*_gvcf.gz"), path("*_gvcf.gz.tbi")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" HaplotypeCaller \
    -R $genome_fasta \
    -I ${bam} \
    -O ${bam.simpleName}_${window.replace(':','_')}_gvcf.gz \
    -ERC GVCF \
    ${gatk_HaplotypeCaller_params}

  """

  stub:
  """
  touch ${bam.simpleName}_gvcf.gz
  """
}

process CombineGVCFs {
  tag "ConsolidateGVCFs"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple path(gvcf), path(gvcf_tbi), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  tuple path("cohort.g.vcf.gz"), path("cohort.g.vcf.gz.tbi")

  script:
  """
  #! /usr/bin/env bash
  GVCFFILES=`echo $gvcf | sed 's/ / --variant /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options "${java_options}" CombineGVCFs \
    -R $genome_fasta \
    --variant \$GVCFFILES \
    --output cohort.g.vcf.gz
  """

  stub:
  """
  #! /usr/bin/env bash
  touch cohort.g.vcf.gz
  """
}

// Joint-Call Cohorts
process GenotypeGVCFs {
  tag "JointCallCohorts"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple path(all_combined_gvcf), path(all_combined_gvcf_tbi), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  tuple path("output.vcf.gz"), path("output.vcf.gz.tbi")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" GenotypeGVCFs \
    -R $genome_fasta \
    -V $all_combined_gvcf \
    --output output.vcf.gz
  """

  stub:
  """
  #! /usr/bin/env bash
  touch output.vcf.gz
  """
}

process calc_DPvalue {
  tag "${vcf_gz.fileName}"
  label 'datamash'

  input:  // sorted SNP vcf
  path(vcf_gz)

  output: // DP value (number) to filter vcf in downstream process
  stdout()

  script:
  """
  #! /usr/bin/env bash
  zcat $vcf_gz | grep -v "^#" | cut -f 8 | grep -oe ";DP=.*" | cut -f 2 -d ';' | cut -f 2 -d "=" > dp.txt
  cat dp.txt | $datamash_app mean 1 sstdev 1 > dp.stats
  cat dp.stats | awk '{print \$1+5*\$2}'
  """

  stub:
  """
  #! /usr/bin/env bash
  # idk what default number should be here
  echo "10"
  """
}

process VariantFiltration {
  tag "${vcf_gz.fileName}"
  publishDir "$params.outdir/05_FilterSNPs", mode: 'copy'

  input:  // [sorted snp vcf, DP filter, genome files ... ]
  tuple path(vcf_gz), path(vcf_gz_tbi), val(dp), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // filtered to identified SNP variants
  path("${vcf_gz.baseName}.marked.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" VariantFiltration \
    --reference $genome_fasta \
    --sequence-dictionary $genome_dict \
    --variant $vcf_gz \
    --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > $dp\" \
    --filter-name "FAIL" \
    --output ${vcf_gz.baseName}.marked.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${vcf_gz.baseName}.marked.vcf
  """
}

