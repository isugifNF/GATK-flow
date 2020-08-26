#!/bin/bash
# script to generate GATK commands for snp calling
# gatk haplotype caller
# for each windows specified
# Arun Seetharam
# 5/16/2019

if [ $# -lt 3 ] ; then
   echo -e "usage: $(basename "$0") <windows bed file> <genome fasta file> <regex for bam files>"
   echo ""
   echo -e "\twindows.bed:\tBED formatted reference intervals to call SNPs"
   echo -e "\tgenome.fasta:\tFASTA formatted reference genome"
   echo -e "\tfinal.bam files:\tBAM formatted final files from processing step"
   echo ""
   exit 0
fi

unset -v bamfiles
list="$1"
bamfiles=(*"${3}")
REF="$3"
for bam in ${bamfiles[@]}; do
echo -en "-I ${bam} ";
done > CombinedBAM_temp


while read line; do
if ! [[ $line == @* ]]; then
g2=$(echo $line | awk '{print $1":"$2"-"$3}'); \
g1=$(echo $line | awk '{print $1"_"$2"_"$3}'); \
CWD=$(pwd)
echo "gatk --java-options \"-Xmx80g -XX:+UseParallelGC\" HaplotypeCaller -R ${REF} $(cat CombinedBAM_temp) -L "${g2}" --output "${g1}".vcf;";
fi
done<${list}
