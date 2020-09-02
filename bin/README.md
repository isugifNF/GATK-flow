# Local instructions of running pipeline

These are temporary notes...for the bash scripts. After the modules are polished then we can delete these notes.

Pull input files at one level up in the `test-data`. This assumes we're running on Condo HPC. 

### Input:

```
test-data/
 |_ ref/     # One reference file
 |  |_ b73_chr1_150000001-151000000.fasta ahalleri
 |
 |_ fastq/   # Several fastq files
    |_ BioSample*_R1.fastq.gz
    |_ BioSample*_R2.fastq.gz
```

### Condo HPC Commands

```
module load samtools
module load picard
module load bwa
module load bedtools2
module load bioawk
bash bin/gatk0_index.sh test-data/ref/b73_chr1_150000001-151000000.fasta ahalleri

module load picard
module load bwa
module load samtools
bash bin/gatk1_preprocess.sh ahalleri.fasta BioSample01 test-data/fastq/BioSample01_R1.fastq.gz test-data/fastq/BioSample01_R2.fastq.gz

bash bin/gatk2b_cmdsgen.sh ahalleri_coords.bed ahalleri.fasta *.bam > gatk3_script.sh

module load gatk
bash gatk3_script.sh
# Similar to: gatk --java-options "-Xmx80g -XX:+UseParallelGC" HaplotypeCaller -R ahalleri.fasta -I BioSample01_R1_final.bam  -L chr1:1-999999 --output chr1_1-999999.vcf
```

The last filter step relies on datamash, which can be downloaded and installed locally ([link to GNU datamash download instructions](https://www.gnu.org/software/datamash/download/)). 

```
module load vcftools
module load gatk
module load bcftools
module load paralleli
bash bin/gatk4_filter_mod.sh ahalleri.fasta

```
