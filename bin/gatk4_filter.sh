#!/bin/bash
# script to filter snps
# gatk tutorial
# Arun Seetharam
# 5/16/2019

merged=merged
#change name as you wish
ref=$1
if [ $# -lt 1 ] ; then
   echo -e "usage: $(basename "$0")  <genome fasta file>"
   echo ""
   echo -e "\tgenome.fasta:\tFASTA formatted reference genome"
   echo ""
   exit 0
fi
module load vcftools
module load GIF/datamash
module load gatk
module load bcftools
mkdir -p vcffiles idxfiles
# merge vcf files
mv *.vcf ./vcffiles
mv *.idx ./idxfiles
cd vcffiles
vcf=(*.vcf)
module load parallel
parallel "grep -v '^#' {}" ::: *.vcf >> ../${merged}.body
grep "^#" ${vcf[1]} > ../${merged}.head
cd ..
cat ${merged}.head ${merged}.body >> ${merged}.vcf
cat ${merged}.vcf | vcf-sort -t $TMPDIR -p 36 -c > ${merged}_sorted.vcf
# calculate stats
bcftools stats ${merged}_sorted.vcf > ${merged}_sorted.vchk
plot-vcfstats ${merged}_sorted.vchk -p plots/
maxdepth=$(grep -oh ";DP=.*;" ${merged}_sorted.vcf | cut -d ";" -f 2 | cut -d "="  -f 2 | datamash mean 1 sstdev 1 | awk '{print $1+$2*5}')
# separate SNPs and INDELs
vcftools --vcf ${merged}_sorted.vcf --keep-only-indels --recode --recode-INFO-all --out ${merged}_sorted-indels.vcf
vcftools --vcf ${merged}_sorted.vcf --remove-indels --recode --recode-INFO-all --out ${merged}_sorted-snps.vcf
gatk --java-options \"-Xmx80g -XX:+UseParallelGC\" VariantFiltration \
     --reference $ref \
    --variant ${merged}_sorted-snps.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 45.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > ${maxdepth}" \
    --filter-name "FAIL" \
    --output ${merged}_filtered-snps.vcf
gatk --java-options \"-Xmx80g -XX:+UseParallelGC\" VariantFiltration \
    --reference $ref \
    --variant ${merged}_sorted-indels.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "FAIL" \
    --output ${merged}_filtered-indels.vcf
