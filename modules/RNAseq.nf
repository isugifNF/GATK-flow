#! /usr/bin/env nextflow
nextflow.enable.dsl=2

process SamToFastq {
  tag "${bam.fileName}"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // reads.bam
  path(bam)

  output: // reads_interleaved.fq
  tuple val("${bam.simpleName}"), path("${bam.simpleName}_newR1.fq.gz"), path("${bam.simpleName}_newR2.fq.gz")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" SamToFastq \
    --INPUT $bam \
    --FASTQ ${bam.simpleName}_newR1.fq.gz \
    --SECOND_END_FASTQ ${bam.simpleName}_newR2.fq.gz \
    --VALIDATION_STRINGENCY SILENT \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_newR1.fq.gz
  touch ${bam.simpleName}_newR2.fq.gz
  """
}

process STAR_index {
  tag "${genome_fasta.simpleName}"
  label 'star'
  publishDir "${params.outdir}/02_MapReads"
  input: tuple path(genome_fasta), path(gtf)

  output: // [genome.fasta, [genome_index files]]
  tuple path("$genome_fasta"), path("Genome/STAR_index")

  script:
  """
  #! /usr/bin/env bash
  $star_app \
  --runThreadN $task.cpus \
  --runMode genomeGenerate \
  --genomeDir Genome/STAR_index \
  --genomeFastaFiles $genome_fasta \
  --sjdbGTFfile $gtf \
  $star_index_params
  """
}

process STAR_align {
  tag "${readname}"
  label 'star'
  publishDir "${params.outdir}/02_MapReads"
  input:
  tuple path(genome_fasta), path(genome_index), val(readname), path(readpairs)

  output:
  tuple val("$readname"), path("star_twopass_output/${readname}_*.bam"), path("star_twopass_output/${readname}_*final.out") // bam? bai?

  script:
  """
  #! /usr/bin/env bash
  mkdir star_twopass_output
  $star_app \
  --runThreadN $task.cpus \
  --outFileNamePrefix star_twopass_output/${readname}_ \
  --genomeDir ${genome_index} \
  --readFilesIn $readpairs \
  --readFilesCommand "gunzip -c" \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic

  # Another option
  # https://gatk.broadinstitute.org/hc/en-us/community/posts/15104189520283-STAR-and-GATK-RNAseq-based-SNP-detection
  """
}

process MergeBamAlignment {
  tag "$i_readname"
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"

  input:  // [readgroup, unmapped reads, mapped reads]
  tuple val(i_readname), path(read_unmapped), path(read_mapped), path(genome_fasta), path(genome_dict)

  output: // merged bam and bai files
  tuple path("${i_readname}_merged.bam"), path("${i_readname}_merged.bai")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" MergeBamAlignment \
  --REFERENCE_SEQUENCE $genome_fasta \
  --UNMAPPED_BAM ${read_unmapped} \
  --ALIGNED_BAM ${read_mapped} \
  --OUTPUT ${i_readname}_merged.bam \
  --INCLUDE_SECONDARY_ALIGNMENTS true \
  --VALIDATION_STRINGENCY SILENT \
  --USE_JDK_DEFLATER true \
  --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${i_readname}_merged.bam
  touch ${i_readname}_merged.bai
  """
}