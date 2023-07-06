#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process SamToFastq {
  tag "${bam.fileName}"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // reads.bam
  path(bam)

  output: // reads_interleaved.fq
  tuple val("${bam.simpleName}"), path("${bam.simpleName}_newR1.fq"), path("${bam.simpleName}_newR2.fq")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" SamToFastq \
    --INPUT $bam \
    --FASTQ ${bam.simpleName}_newR1.fq \
    --SECOND_END_FASTQ ${bam.simpleName}_newR2.fq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INCLUDE_NON_PF_READS true \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_newR1.fq
  touch ${bam.simpleName}_newR2.fq
  """
}

// INTERLEAVE=true
// USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

process bwamem2_index {
  tag "${genome_fasta.simpleName}"
  label 'bwamem'
  publishDir "${params.outdir}/02_MapReads"

  input:
  path(genome_fasta)

  output: // [genome.fasta, [genome_index files]]
  tuple path("$genome_fasta"), path("${genome_fasta}*")

  script:
  """
  #! /usr/bin/env bash
  $bwamem2_app index $genome_fasta
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta}.variousindexfiles
  """
}

process bwamem2_mem {
  tag "$readname"
  label 'bwamem'
  publishDir "${params.outdir}/02_MapReads"

  input:
  tuple path(genome_fasta), path(genome_index), val(readname), path(readpairs)

  output: // reads_mapped_2_genome.bam
  tuple val("$readname"), path("${readname}_mapped.bam")

  script:
  """
  #! /usr/bin/env bash
  PROC1=\$((`nproc` * 3/4))
  $bwamem2_app mem -t \${PROC1} ${genome_fasta} ${readpairs} |\
     $samtools_app view --threads 1 -bS - > ${readname}_mapped.bam
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${readname}_mapped.bam
  """
}