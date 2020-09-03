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
```

<details><summary>See list of generated files</summary>

```
-rw-r--r--. 1 jenchang its-hpc-condo-severin 1012505 Sep  3 11:23 ahalleri.fasta
-rw-r--r--. 1 jenchang its-hpc-condo-severin      12 Sep  3 11:23 ahalleri.length
-rw-r--r--. 1 jenchang its-hpc-condo-severin  500048 Sep  3 11:23 ahalleri.fasta.sa
-rw-r--r--. 1 jenchang its-hpc-condo-severin  250001 Sep  3 11:23 ahalleri.fasta.pac
-rw-r--r--. 1 jenchang its-hpc-condo-severin      20 Sep  3 11:23 ahalleri.fasta.fai
-rw-r--r--. 1 jenchang its-hpc-condo-severin 1000072 Sep  3 11:23 ahalleri.fasta.bwt
-rw-r--r--. 1 jenchang its-hpc-condo-severin      37 Sep  3 11:23 ahalleri.fasta.ann
-rw-r--r--. 1 jenchang its-hpc-condo-severin      11 Sep  3 11:23 ahalleri.fasta.amb
-rw-r--r--. 1 jenchang its-hpc-condo-severin     146 Sep  3 11:23 ahalleri.dict
-rw-r--r--. 1 jenchang its-hpc-condo-severin      14 Sep  3 11:23 ahalleri_coords.list
```
 
</details>

```
module load picard
module load bwa
module load samtools
bash bin/gatk1_preprocess.sh ahalleri.fasta BioSample01 test-data/fastq/BioSample01_R1.fastq.gz test-data/fastq/BioSample01_R2.fastq.gz
```

<details><summary>See list of generated files</summary>
 
```
-rw-r--r--. 1 jenchang its-hpc-condo-severin      946 Sep  3 11:29 BioSample01_R1_markilluminaadapters_metrics.txt
-rw-r--r--. 1 jenchang its-hpc-condo-severin 92652492 Sep  3 11:30 BioSample01_R1_prefinal.bam
-rw-r--r--. 1 jenchang its-hpc-condo-severin     3304 Sep  3 11:30 BioSample01_R1_prefinal.bai
-rw-r--r--. 1 jenchang its-hpc-condo-severin     2959 Sep  3 11:30 BioSample01_R1_mergebamalignment_markduplicates_metrics.txt
-rw-r--r--. 1 jenchang its-hpc-condo-severin 92576481 Sep  3 11:30 BioSample01_R1_final.bam
-rw-r--r--. 1 jenchang its-hpc-condo-severin     3352 Sep  3 11:30 BioSample01_R1_final.bai
```

</details>

```
bash bin/gatk2b_cmdsgen.sh ahalleri_coords.list ahalleri.fasta *_final.bam > gatk3_script.sh
```

Which generates gatk commands and puts them in a `gatk3_script.sh` script.

<details><summary>See example gatk3_script.sh</summary>
 
**gatk3_script.sh**

```
gatk --java-options "-Xmx80g -XX:+UseParallelGC" HaplotypeCaller -R ahalleri.fasta -I BioSample01_R1_final.bam  -L chr1:1-999999 --output chr1_1-999999.vcf;
```
 
</details>

```
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
