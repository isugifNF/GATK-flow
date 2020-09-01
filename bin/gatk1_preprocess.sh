#!/bin/bash
# script to prepare fastq files for GATK snp calling
# Arun Seetharam
# 5/16/2019
if [ $# -ne 3 ] ; then
   echo -e "usage: $(basename "$0") <reference> <read-group> <read R1> <read R2>"
   echo ""
   echo -e "\treference:\tindexed reference genome (full path)-short name in indexing script is fine"
   echo -e "\tread-group:\t the line name or individual name for which the fastq files belong to"
   echo -e "\tR1:\t forward read"
   echo -e "\tR2:\t reverse read"
   echo ""
   exit 0
fi
module load picard
module load bwa
module load samtools
ulimit -c unlimited
REF=$1
RGID=$2
R1=$3
R2=$4
# adjust this to suit your input file name
OUT=$(basename ${R1%%.*}) # stripping fastq.gz from R1 file
PICARD_HOME=$(dirname $(which picard))
PICARD_CMD="java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_HOME}/picard.jar"

# platform id from fastq file
if [ ${R1: -3} == ".gz" ]; then
   PLT=$(zcat $R1 |head -n 1 |cut -f 3 -d ":")
   RGPU=$(zcat $R1 |head -n 1 |cut -f 3-5 -d ":")
else
   PLT=$(cat $R1 |head -n 1 |cut -f 3 -d ":")
   RGPU=$(cat $R1 |head -n 1 |cut -f 3-5 -d ":")
fi

# time stamp as string of numbers
TDATE=$(date '+%Y-%m-%d %H:%M:%S' |sed 's/ /T/g')
#  library identifier
RGLB="$RGID"
# platform name choose either ILLUMINA, SOLID, LS454, HELICOS and PACBIO
RGPL="ILLUMINA"
# genotype name, this will appear in VCF file header
RGSM="$RGID"
# convert fastq to sam and add readgroups
$PICARD_CMD FastqToSam \
   FASTQ=${R1} \
   FASTQ2=${R2} \
   OUTPUT=${OUT}_fastqtosam.bam \
   READ_GROUP_NAME=${OUT} \
   SAMPLE_NAME=${OUT}_name \
   LIBRARY_NAME=${OUT}_lib \
   PLATFORM_UNIT=${PLT} \
   PLATFORM=illumina \
   SEQUENCING_CENTER=ISU \
   RUN_DATE=${TDATE}  || {
echo >&2 ERROR: FastqToSam failed for $OUT
exit 1
}
# marking adapters
$PICARD_CMD MarkIlluminaAdapters \
   I=${OUT}_fastqtosam.bam \
   O=${OUT}_markilluminaadapters.bam \
   M=${OUT}_markilluminaadapters_metrics.txt \
   TMP_DIR=${TMPDIR}  || {
echo >&2 ERROR: MarkIlluminaAdapters failed for $OUT
exit 1
}
# convert bam back to fastq for mapping
$PICARD_CMD SamToFastq \
   I=${OUT}_markilluminaadapters.bam \
   FASTQ=${OUT}_samtofastq_interleaved.fq \
   CLIPPING_ATTRIBUTE=XT \
   CLIPPING_ACTION=2 \
   INTERLEAVE=true NON_PF=true TMP_DIR=${TMPDIR} || {
echo >&2 ERROR: SamToFastq failed for $OUT
exit 1
}
# mapping reads to indexed genome
bwa mem \
   -M \
   -t 15 \
   -p $REF \
   ${OUT}_samtofastq_interleaved.fq |\
  samtools view -buS - > ${OUT}_bwa_mem.bam || {
echo >&2 ERROR: BWA failed for $OUT
exit 1
}
# merging alignments
$PICARD_CMD MergeBamAlignment \
   R=$REF \
   UNMAPPED_BAM=${OUT}_fastqtosam.bam \
   ALIGNED_BAM=${OUT}_bwa_mem.bam \
   O=${OUT}_snippet_mergebamalignment.bam \
   CREATE_INDEX=true \
   ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \
   CLIP_OVERLAPPING_READS=true \
   INCLUDE_SECONDARY_ALIGNMENTS=true \
   MAX_INSERTIONS_OR_DELETIONS=-1 \
   PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
   ATTRIBUTES_TO_RETAIN=XS \
   TMP_DIR="${TMPDIR}" || {
echo >&2 ERROR: MergeBamAlignment failed for $OUT
exit 1
}
# mark duplicates
$PICARD_CMD MarkDuplicates \
  INPUT=${OUT}_snippet_mergebamalignment.bam \
  OUTPUT=${OUT}_prefinal.bam \
  METRICS_FILE=${OUT}_mergebamalignment_markduplicates_metrics.txt \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  CREATE_INDEX=true \
  TMP_DIR=$TMPDIR || {
echo >&2 ERROR: MarkDuplicates failed for $OUT
exit 1
}
# add read groups
$PICARD_CMD AddOrReplaceReadGroups \
  INPUT=${OUT}_prefinal.bam \
  OUTPUT=${OUT}_final.bam \
  RGID=${RGID} \
  RGLB=${RGLB} \
  RGPL=${RGPL} \
  RGPU=${RGPU} \
  RGSM=${RGSM} \
  CREATE_INDEX=true \
  TMP_DIR=$TMPDIR || {
echo >&2 ERROR: Adding read groups failed for $OUT
exit 1
}
echo >&2 "ALL DONE!"
# cleanup
rm ${OUT}_fastqtosam.bam
rm ${OUT}_markilluminaadapters.bam
rm ${OUT}_samtofastq_interleaved.fq
rm ${OUT}_bwa_mem.bam
rm ${OUT}_snippet_mergebamalignment.bam
rm ${OUT}_snippet_mergebamalignment.bai
