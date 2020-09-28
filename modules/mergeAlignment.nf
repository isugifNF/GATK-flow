#! /usr/bin/env nextflow

picard_container = ''

process MergeBamAlignment_help {
  label 'picard'

  container = "$picard_container"

  output: path 'picard-MergeBamAlignment.txt'

  """
  picard MergeBamAlignment > picard-MergeBamAlignment.txt
  """
}

process test_run {
    tag "$readname"
    label 'picard'
    publishDir "${params.outdir}/MergeBamAlignment", mode: 'copy'

    input:
    tuple val(readname), path(read_unmapped), path(read_aligned), path(genome_fasta), path(genome_index), path(genome_coords)

    output:
    stdout()

    script:
    """
    #! /usr/bin/env bash

    echo "readname = $readname"
    echo "read_unmapped = $read_unmapped"
    echo "read_aligned = $read_aligned"
    echo "genome_fasta = $genome_fasta"
    echo "genome_index = $genome_index"
    echo "genome_coords = $genome_coords"
    """
}

process MergeBamAlignment_run {
    tag "$readname"
    label 'picard'
    publishDir "${params.outdir}/MergeBamAlignment", mode: 'copy'

    input:
    tuple val(readname), path(read_unmapped), path(read_aligned), path(genome_fasta), path(genome_index), path(genome_coords), path(fai)

    output:
    path "${readname}_merged.bam"
    path "${readname}_merged.bai"

    script:
    """
    #! /usr/bin/env bash

    picard MergeBamAlignment \
    R=$genome_fasta \
    UNMAPPED_BAM=$read_unmapped \
    ALIGNED_BAM=$read_aligned \
    O=${readname}_merged.bam \
    CREATE_INDEX=true \
    ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false \
    CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    ATTRIBUTES_TO_RETAIN=XS
    """
}

/* Scrap here
picard MergeBamAlignment \
R=$genome_fasta \
UNMAPPED_BAM=${read_unmapped} \
ALIGNED_BAM=${read_aligned} \
O=${readname}_MergeBamAlignment.bam \
CREATE_INDEX=true \
ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false \
CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true \
MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
ATTRIBUTES_TO_RETAIN=XS
*/
