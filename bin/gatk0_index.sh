#!/bin/bash
# script to prepare reference genome for GATK snp calling
# indexes for BWA mapping porgram
# generates windows
# Arun Seetharam
# 5/16/2019

if [ $# -ne 2 ] ; then
   echo -e "usage: $(basename "$0") <genome.fasta> <name>"
   echo ""
   echo -e "\tgenome.fasta:\tFASTA formatted reference genome"
   echo -e "\tname:\tsmall name for easy reference management, any argument will suffice, must be one word"
   echo ""
   exit 0
fi

module load samtools
module load picard
module load bwa
module load bedtools2
module load bioawk
ref="$1"
name="$2"
window=10000000
bioawk -c fastx '{print}' $ref | sort -k1,1V | awk '{print ">"$1;print $2}' | fold > ${name}.fasta
picard CreateSequenceDictionary REFERENCE=${name}.fasta OUTPUT=${name}.dict
bwa index -a bwtsw ${name}.fasta
bioawk -c fastx '{print $name"\t"length($seq)}' ${name}.fasta > ${name}.length
bedtools makewindows -w $window -g ${name}.length |\
   awk '{print $1"\t"$2+1"\t"$3}' |\
   sed 's/\t/:/1' |\
   sed 's/\t/-/1' > ${name}_coords.bed
