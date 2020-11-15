#! /usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMsg() {
  log.info """
   Usage:
   The typical command for running the pipeline is as follows:
   nextflow run main.nf --genome GENOME.fasta --reads "*_{R1,R2}.fastq.gz" -profile singularity
   nextflow run main.nf --genome GENOME.fasta --reads_file READ_PATHS.txt -profile singularity

   Mandatory arguments:
    --genome                Genome fasta file, against which reads will be mapped to find SNPs
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
    or
    --genome                Genome fasta file, against which reads will be mapped to find SNPs
    --reads_file            Text file (tab delimited) with three columns [readname left_fastq.gz right_fastq.gz]. Will need full path for files.

   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, condo, atlas, singularity [default:local]
    --singularity_img       Singularity image if [-profile singularity] is set [default:'shub://aseetharam/gatk:latest']
    --bwa_app               Link to bwa executable [default: 'bwa']
    --samtools_app          Link to samtools executable [default: 'samtools']
    --picard_app            Link to picard executable [default: 'picard'], might want to change to "java -jar ~/PICARD_HOME/picard.jar"
    --bedtools_app          Link to bedtools executable [default: 'bedtools']
    --gatk_app              Link to gatk executable [default: 'gatk']
    --datamash_app          Link to datamash executable [default: 'datamash']
    --vcftools_app          Link to vcftools executable [default: 'vcftools']

   Optional other arguments:
    --window                Window size passed to bedtools for gatk [default:100000]
    --queueSize             Maximum jobs to submit to slurm [default:18]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help
    
"""
}

if(params.help){
  helpMsg()
  exit 0
}

if(!params.genome) {
  log.info"""
#===============
  ERROR: --genome GENOME.fasta    A reference genome file is required!
#===============
  """
  helpMsg()
  exit 0
}

if(!params.reads & !params.reads_file){
  log.info"""
#===============
  ERROR: --reads "*_{r1,r2}.fq.gz"     Paired-end read files are required! Either as a glob or as a tab-delimited text file
         --reads_file READS_FILE.txt
#===============
  """
  helpMsg()
  exit 0
}

process get_test_data {
  publishDir './', mode: 'move'

  output:
  path "*"

  script:
  """
  #! /usr/bin/env bash
  wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
  tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
  """
}

process fasta_sort {
  tag "$fasta"
  publishDir "$params.outdir/sort_fasta"


  input:
  path fasta

  output:
  path "${fasta.simpleName}_sorted.fasta"

  script:
  """
  #! /usr/bin/env bash
  cat $fasta |\
    tr '\n' '\t' | sed \$'s/>/\r>/g'| tr '\r' '\n'|\
    sort |\
    tr '\t' '\n' | sed \$'s/\r//g' |\
    grep -v "^\$" > ${fasta.simpleName}_sorted.fasta
  """
}

process fasta_bwa_index {
  tag "$fasta"
  label 'bwa'
  publishDir "${params.outdir}/bwa"

  input:
  path fasta

  output: // [genome.fasta, [genome_index files]]
  tuple path("$fasta"), path("${fasta}*")

  script:
  """
  #! /usr/bin/env bash
  $bwa_app index $fasta
  """
}

process fasta_samtools_faidx {
  tag "$fasta"
  label 'samtools'
  publishDir "${params.outdir}/samtools"

  input:
  path fasta

  output:
  path "${fasta}.fai"

  """
  #! /usr/bin/env bash
  $samtools_app faidx $fasta
  """
}

process fasta_picard_dict {
  tag "$fasta"
  label 'picard'
  publishDir "${params.outdir}/picard"

  input:
  path fasta

  output:
  path "${fasta.simpleName}.dict"

  script:
  """
  #! /usr/bin/env bash
  $picard_app CreateSequenceDictionary \
    REFERENCE=${fasta} \
    OUTPUT=${fasta.simpleName}.dict
  """
}

process paired_FastqToSAM {
  tag "$readname"
  label 'picard'
  publishDir "${params.outdir}/picard"

  input:
  tuple val(readname), path(readpairs)

  output:
  path "${readname}.bam"

  """
  #! /usr/bin/env bash
  $picard_app FastqToSam \
    FASTQ=${readpairs.get(0)} \
    FASTQ2=${readpairs.get(1)} \
    OUTPUT=${readname}.bam \
    READ_GROUP_NAME=${readname} \
    SAMPLE_NAME=${readname}_name \
    LIBRARY_NAME=${readname}_lib \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="ISU" \
    USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
  """
}

process BAM_MarkIlluminaAdapters {
  tag "${bam.fileName}"
  label 'picard'
  publishDir "${params.outdir}/picard"

  input:
  path bam

  output:
  path "${bam.simpleName}_marked.bam"  //, emit: bam
  //path "${bam.simpleName}_marked*.txt"

  script:
  """
  #! /usr/bin/env bash
  $picard_app MarkIlluminaAdapters \
    I=$bam \
    O=${bam.simpleName}_marked.bam \
    M=${bam.simpleName}_marked_metrics.txt \
    USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
  """
}

process BAM_SamToFastq {
  tag "${bam.fileName}"
  label 'picard'
  publishDir "${params.outdir}/picard"

  input:
  path bam

  output:
  path "${bam.simpleName}_interleaved.fq"

  script:
  """
  #! /usr/bin/env bash
  $picard_app SamToFastq \
    I=$bam \
    FASTQ=${bam.simpleName}_interleaved.fq \
    CLIPPING_ATTRIBUTE=XT \
    CLIPPING_ACTION=2 \
    INTERLEAVE=true \
    NON_PF=true \
    USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
  """
}

process run_bwa_mem {
  tag "$readsfq"
  label 'bwa_mem'
  publishDir "${params.outdir}/bwa_mem"

  input:
  tuple path(readsfq), path(genome_fasta), path(genome_index)

  output:
  path "${readsfq.simpleName}_mapped.bam"

  script:
  """
  #! /usr/bin/env bash
  $bwa_app mem \
   -M \
   -t 15 \
   -p ${genome_fasta} \
   ${readsfq} |\
  samtools view -buS - > ${readsfq.simpleName}_mapped.bam
  """
}

process run_MergeBamAlignment {
  tag "$readname"
  label 'picard'
  publishDir "${params.outdir}/picard"

  input:
  tuple val(readname), path(read_unmapped), path(read_mapped), path(genome_fasta), path(genome_index), path(genome_fai), path(genome_dict)

  output:
  path "${readname}_merged.bam", emit:'bam'
  path "${readname}_merged.bai", emit:'bai'

  script:
  """
  #! /usr/bin/env bash
  $picard_app MergeBamAlignment \
    R=$genome_fasta \
    UNMAPPED_BAM=$read_unmapped \
    ALIGNED_BAM=$read_mapped \
    O=${readname}_merged.bam \
    CREATE_INDEX=true \
    ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false \
    CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    ATTRIBUTES_TO_RETAIN=XS \
    USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
  """
}

params.window=100000

process fai_bedtools_makewindows {
  tag "$fai"
  label 'bedtools'
  publishDir "${params.outdir}/bedtools"

  input:
  path fai

  output:
  path "${fai.simpleName}_coords.bed"

  script:
  """
  #! /usr/bin/env bash
  awk -F'\t' '{print \$1"\t"\$2}' $fai > genome_length.txt
  $bedtools_app makewindows -w $params.window -g genome_length.txt |\
    awk '{print \$1"\t"\$2+1"\t"\$3}' |\
    sed \$'s/\t/:/1' |\
    sed \$'s/\t/-/g' > ${fai.simpleName}_coords.bed
  """
}

process run_gatk_snp {
  tag "$window"
  label 'gatk'
  publishDir "${params.outdir}/gatk"

  input:
  tuple val(window), path(bam), path(bai), path(genome_fasta), path(genome_index), path(genome_fai), path(genome_dict)

  output:
  path "*.vcf", emit: 'vcf'
  path "*.vcf.idx", emit: 'idx'

  script:
  """
  #! /usr/bin/env bash
  BAMFILES=`echo $bam | sed 's/ / -I /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options \"-Xmx80g -XX:+UseParallelGC\" HaplotypeCaller -R $genome_fasta -I \$BAMFILES -L $window --output ${window.replace(':','_')}.vcf
  """
}

process merge_vcf {
  publishDir "${params.outdir}/vcftools"

  input:
  path(vcfs)

  output:
  path "first-round_merged.vcf"

  script:
  """
  #! /usr/bin/env bash
  cat ${vcfs.get(0)} | grep "^#" > first-round_merged.vcf
  cat ${vcfs} | grep -v "^#" >> first-round_merged.vcf
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
  tag "$vcf.fileName"
  label 'picard'
  publishDir "$params.outdir/picard"

  input:
  tuple path(vcf), path(dict)

  output:
  tuple path ("*.vcf"), path("*.vcf.*")

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
  cat dp.txt | $datamash_app mean 1 sstdev 1 > dp.stats
  cat dp.stats | awk '{print \$1+5*\$2}'
  """
}

process gatk_VariantFiltration {
  tag "$sorted_snp_vcf.fileName"

  input:
  tuple path(sorted_snp_vcf), val(dp), path(genome_fasta), path(genome_index), path(genome_fai), path(genome_dict)

  output:
  path("${sorted_snp_vcf.simpleName}.marked.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app VariantFiltration \
    --reference $genome_fasta \
    --sequence-dictionary $genome_dict \
    --variant $sorted_snp_vcf \
    --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > $dp\" \
    --filter-name "FAIL" \
    --output ${sorted_snp_vcf.simpleName}.marked.vcf
  """
}

// java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar

process keep_only_pass {
  tag "$snp_marked_vcf.fileName"

  input:
  path(snp_marked_vcf)

  output:
  path("${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf")

  script:
  """
  #! /usr/bin/env bash
  cat $snp_marked_vcf | grep "^#" > ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  cat $snp_marked_vcf | grep -v "^#" | awk '\$7=="PASS"' >> ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  """
}

workflow prep_genome {
  take: reference_fasta
  main:
    fasta_sort(reference_fasta) | (fasta_bwa_index & fasta_samtools_faidx & fasta_picard_dict )

    genome_ch = fasta_bwa_index.out
      .combine(fasta_samtools_faidx.out)
      .combine(fasta_picard_dict.out)

  emit:
    genome_ch
}

workflow prep_reads {
  take: reads_fastas
  main:
    reads_fastas | paired_FastqToSAM | BAM_MarkIlluminaAdapters //| BAM_SamToFastq

    reads_ch = BAM_MarkIlluminaAdapters.out

  emit:
    reads_ch
}

workflow map_reads {
   take:
     reads_ch
     genome_ch

   main:
     bwaindex_ch = genome_ch.map { n -> [ n.get(0), n.get(1) ] }
     reads_ch | BAM_SamToFastq
     BAM_SamToFastq.out.combine(bwaindex_ch) | run_bwa_mem

     reads_mapped_ch = run_bwa_mem.out.map { n -> [ n.simpleName.replaceFirst("_marked_interleaved_mapped", ""), n ] }
     reads_unmapped_ch = reads_ch.map { n -> [ n.simpleName.replaceFirst("_marked",""), n ] }

     reads_merged_ch = reads_unmapped_ch
       .join(reads_mapped_ch)
       .combine(genome_ch)

   emit:
     reads_merged_ch
}

workflow {
  main:
    if ("$workflow.profile"=~/testdata/) {
      get_test_data()
    } else {
      genome_ch = channel.fromPath(params.genome, checkIfExists:true) | prep_genome // | view
    if (params.reads) {
      reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true) | prep_reads //| view
    } else {
      reads_ch = channel.fromPath(params.reads_file, checkIfExists:true).splitCsv(sep:'\t') |
         map { n -> [ n.get(0), [ n.get(1), n.get(2) ]] } | prep_reads
    }

    map_reads(reads_ch, genome_ch) | run_MergeBamAlignment

    bam_ch = run_MergeBamAlignment.out.bam.toList() | map { it -> [it]} //| view
    bai_ch = run_MergeBamAlignment.out.bai.toList() | map { it -> [it]} //| view
    merged_bam_ch = bam_ch.combine(bai_ch) //| view

    genome_ch.map { n -> n.get(2) } |
      fai_bedtools_makewindows |
      splitText( by:1 ) |
      map { n -> n.replaceFirst("\n","") } |
      combine(merged_bam_ch) |
      combine(genome_ch) |
      run_gatk_snp

    run_gatk_snp.out.vcf.toList() |
      merge_vcf |
      vcftools_snp_only |
      combine(genome_ch.map {n -> n.get(3)}) |
      run_SortVCF |
      map{n->n.get(0)} |
      calc_DPvalue |
      view

    run_SortVCF.out.map{n -> n.get(0)} |
      combine(calc_DPvalue.out.map{n-> n.replaceAll("\n","")}) |
      combine(genome_ch) |
      gatk_VariantFiltration |
      keep_only_pass |
      view
    }

}
