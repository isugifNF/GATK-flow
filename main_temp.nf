#! /usr/bin/env nextflow

nextflow.enable.dsl=2

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

  output:
  path "$fasta", emit: genome_fasta
  path "${fasta}*", emit: genome_index

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

//picard_app='java -jar /picard/picard.jar'
//picard_app='java -jar ~/bin/picard.jar'

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

workflow prep_genome {
  take: reference_fasta
  main:
    fasta_sort(reference_fasta) | (fasta_bwa_index & fasta_samtools_faidx & fasta_picard_dict )

    genome_ch = fasta_bwa_index.out.genome_fasta
      .combine(fasta_bwa_index.out.genome_index.toList())
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
    //get_test_data()
    genome_ch = channel.fromPath(params.genome, checkIfExists:true) | prep_genome // | view
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true).take(3)| prep_reads //| view

    map_reads(reads_ch, genome_ch) | run_MergeBamAlignment

    bam_ch = run_MergeBamAlignment.out.bam.toList() | map { it -> [it]} //| view
    bai_ch = run_MergeBamAlignment.out.bai.toList() | map { it -> [it]} //| view
    merged_bam_ch = bam_ch.combine(bai_ch) //| view


    genome_ch.map { n -> n.get(2) } | fai_bedtools_makewindows
    gatk_input_ch = fai_bedtools_makewindows
      .out.splitText( by:1 )
      .map { n -> n.replaceFirst("\n","") }
      .combine(merged_bam_ch)
      .combine(genome_ch) | run_gatk_snp
}
