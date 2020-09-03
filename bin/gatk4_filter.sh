#!/bin/bash
# filter the SNPs for running BSQR from the first round GATK
# needs the directory of VCF files created from the first round GATK
# Arun Seetharam
# arnstrm@iastate.edu
# 08/23/2020
if [ $# -lt 2 ] ; then
        echo "usage: ${0##*/}  /path/containing/vcf/splitfiles/ genome.fasta"
        echo ""
        echo "Filters the SNPs for running BSQR from the first round GATK"
        echo ""
        echo "1. Full path for the location of VCF files generated from the first round"
		echo "   These files will be combined first and then filtered"
		echo "2. Reference genome used for SNP calling"
		echo "   This assumes that dict file with dict extension is located in the same dir"
        echo ""
        echo ""
		echo "expects modules 'gatk', 'picard', 'vcftools', and executable 'datamash'"
        echo ""
exit 0
fi


vcfdir=$1
ref=$2
dict=${ref%.*}.dict
if [[ ! -f ${dict} ]]
then
    echo "dict file does not exist for the reference file"
fi
# check if datamash is available, if using container, you can probably omit this
cmd=datamash
type -P "$cmd" && echo "$cmd found" || echo "$cmd not found"
# set the dir (tmpdir has fast I/O so we will use that)
workdir=${TMPDIR}
destdir=${SLURM_SUBMIT_DIR}
# concatenate files to create single VCF file from interval VCF files
cd ${workdir}
for vcf in ${vcfdir}/*.vcf; do
	base=$(basename ${vcf%.*});
	grep -v "^#" ${vcf} > ${workdir}/${base}.body
done
grep "^#" $vcf > head.txt
cat *.body >> full.txt
cat head.txt full.txt > first-round_merged.vcf
rnd1=first-round_merged.vcf
# remove indels from the VCF file
ml purge
ml vcftools
vcftools \
    --vcf ${rnd1} \
	--remove-indels \
	--recode \
	--recode-INFO-all \
	--out ${rnd1%.*}_snps-only
# sort the vcf file using picard
# the default picard uses very less -Xmx value and hence we need to request more
# also creates the index file needed for next step
ml picard
picardhome=$(dirname $(which picard))
java -Xmx100g -Djava.io.tmpdir=$TMPDIR  -jar ${picardhome}/picard.jar SortVcf \
	INPUT=${rnd1%.*}_snps-only.recode.vcf \
	SEQUENCE_DICTIONARY ${dict} \
	CREATE_INDEX=true \
	OUTPUT=${rnd1%.*}_snps-only.sorted.vcf
# calculate DP value for filtering
# gatk bp recommends mean + (5*sd)
grep -v "^#" ${rnd1%.*}_snps-only.sorted.vcf | cut -f 8 | grep -oe ";DP=.*" |cut -f 2 -d ";" | cut -f 2 -d "=" > dp.txt
cat dp.txt | datamash mean 1 sstdev 1 > dp.stats
DP=$(cat dp.stats  |awk '{print $1+5*$2}')
# perform filteration using the expression:
# QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > $DP
# will use GATK with extra -Xmx value
ml purge
ml gatk
gatkpath=$(dirname $(which gatk))
gatk=$(basename $(find $gatkdir -name "*local.jar"))
java -Xmx100g -Djava.io.tmpdir=$TMPDIR  -jar ${gatkpath}/${gatk} VariantFiltration \
	--reference $ref \
	--variant ${rnd1%.*}_snps-only.sorted.vcf \
	--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > $DP" \
	--filterName "FAIL" \
	--output ${rnd1%.*}_snps-only.marked.vcf
# retain only the one with PASS. I think we can generate this directly from the previous gatk step,
# but sometimes it is useful to have both versions
grep -v "^#" ${rnd1%.*}_snps-only.marked.vcf | awk '$7=="PASS"' > body.pass
grep "^#" ${rnd1%.*}_snps-only.marked.vcf > head.pass
cat  head.pass body.pass >> ${rnd1%.*}_snps-only.pass-only.vcf
# rsync files from tmp to destination dir
rsync -avP ${workdir}/ ${destdir}/
