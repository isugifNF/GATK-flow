# Maize_WGS_Build
Code and Data to include in WGS Build

Start with the workflow documented on [Bioinformatic Workbook GATK DNAseq Best Practices](https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gatk-dnaseq-best-practices-workflow.html#gsc.tab=0)

### Test dataset

A simple test dataset is available (here)[/test-data]. This dataset contains a small genome (portion of chr1, B73v5), and Illumina short reads for 26 NAM lines (including B73) and B73Ab10 line (27 lines total).
Only the reads that map to the region of the v5 genome is included, so that this can be tested quickly.
There are examples of multiple files belonging to same NAM line as well as single file per NAM line to make sure both conditions works correctly.
The end VCF file should have exactly 27 individuals (lines) in them.

* Test DataSet is on ISU Box -- [https://iastate.app.box.com/v/gatk-test-data](https://iastate.app.box.com/v/gatk-test-data)

### Container 

Tools required for the workflow are included in the container

#### To pull the image

```
singularity pull --name gatk.sif shub://aseetharam/gatk:latest
```

#### To use the image

```
singularity exec gatk.sif samtools
singularity exec gatk.sif bwa
singularity exec gatk.sif datamash
singularity exec gatk.sif gatk
singularity exec gatk.sif java -jar /picard/picard.jar
singularity exec gatk.sif vcftools
```

### Running the pipeline

<!--
>
> If on a local laptop with nextflow installed:
> 
> ```
> nextflow run HuffordLab/Maize_WGS_Build
> ```
> 
> If on HPCC Condo:
> 
> ```
> module load gcc/7.3.0-xegsmw4 nextflow
> nextflow run HuffordLab/Maize_WGS_Build -profile condo
> ```
-->

On HPCC Condo:

```
git clone https://github.com/HuffordLab/Maize_WGS_Build.git
cd Maize_WGS_Build
# Fetch dataset from ISU box here and place in `test-data` folder

module load gcc/7.3.0-xegsmw4 nextflow
nextflow run main.nf -profile condo
```

On MacOS laptop where dependencies are locally installed:

```
$ nextflow run main_temp.nf --picard_app "java -jar ~/bin/picard.jar"

N E X T F L O W  ~  version 20.07.1
Launching `main_temp.nf` [loving_euclid] - revision: 11fbfd945a
executor >  local (24)
[2c/c03319] process > prep_genome:fasta_sort (b73_chr1_150000001-151000000.fasta)                  [100%] 1 of 1 ✔
[98/947da7] process > prep_genome:fasta_bwa_index (b73_chr1_150000001-151000000_sorted.fasta)      [100%] 1 of 1 ✔
[5f/a3dc4f] process > prep_genome:fasta_samtools_faidx (b73_chr1_150000001-151000000_sorted.fasta) [100%] 1 of 1 ✔
[77/4c9fb5] process > prep_genome:fasta_picard_dict (b73_chr1_150000001-151000000_sorted.fasta)    [100%] 1 of 1 ✔
[2e/c31cba] process > prep_reads:paired_FastqToSAM (BioSample19)                                   [100%] 3 of 3 ✔
[39/b2bab3] process > prep_reads:BAM_MarkIlluminaAdapters (BioSample19.bam)                        [100%] 3 of 3 ✔
[56/8606d8] process > map_reads:BAM_SamToFastq (BioSample19_marked.bam)                            [100%] 3 of 3 ✔
[31/6b6834] process > map_reads:run_bwa_mem (BioSample19_marked_interleaved.fq)                    [100%] 3 of 3 ✔
[4d/06ad1f] process > run_MergeBamAlignment (BioSample19)                                          [100%] 3 of 3 ✔
[fc/e839ea] process > fai_bedtools_makewindows (b73_chr1_150000001-151000000_sorted.fasta.fai)     [100%] 1 of 1 ✔
[2c/8fe3fc] process > run_gatk_snp (chr1:400001-500000)                                            [  0%] 0 of 10
```

<details><summary>See generated <b>results</b> folder</summary>

```
ls -l results/
#> total 5736
#> drwxr-xr-x  3 jenchang  staff    96B Sep 10 18:36 bedtools
#> drwxr-xr-x  8 jenchang  staff   256B Sep 10 18:36 bwa
#> drwxr-xr-x  3 jenchang  staff    96B Sep 10 18:36 createSeqDict
#> drwxr-xr-x  3 jenchang  staff    96B Sep 10 18:36 faidx
#> -rw-r--r--  1 jenchang  staff   2.8M Sep 10 18:36 report.html
#> drwxr-xr-x  3 jenchang  staff    96B Sep 10 18:36 seqLength
#> drwxr-xr-x  3 jenchang  staff    96B Sep 10 18:36 sortSeq
#> -rw-r--r--  1 jenchang  staff   6.4K Sep 10 18:36 timeline.html
```

</details>

<!--

<details><summary>See example HPCC Condo running output </summary>

In this case there are 101 slurm jobs on the queue so far. The process `fastqc` has a total of 258 jobs to submit (one for each `test-data` fastq file).

```
nextflow run main.nf -profile condo
#> N E X T F L O W  ~  version 20.07.1
#> Launching `main.nf` [boring_carson] - revision: 99983aad6a
#> executor >  slurm (101)
#> [0f/70feab] process > fastqc (null)         [  0%] 1 of 258
#> [f4/0b666a] process > gatk0_index_help      [  0%] 0 of 1
#> [ef/d2fbd1] process > gatk0_index (1)       [  0%] 0 of 1
#> [2d/c71570] process > gatk2_preprocess_help [100%] 1 of 1 ✔
#> [57/4481cd] process > gatk3_cmdsgen_help    [100%] 1 of 1 ✔
#> [cf/a201a6] process > gatk4_filter_help     [100%] 1 of 1 ✔
#> /work/GIF/jenchang/_wrkspc/_testremote/Maize_WGS_Build/test-data/ref/b73_chr1_150000001-151000000.fasta
#> /work/GIF/jenchang/_wrkspc/_testremote/Maize_WGS_Build/test-data/fastq/1721-5_S1_L004_R1_001.fastq.gz
#> /work/GIF/jenchang/_wrkspc/_testremote/Maize_WGS_Build/test-data/fastq/CML333_S0_L001_R2_001.fastq.gz
#> /work/GIF/jenchang/_wrkspc/_testremote/Maize_WGS_Build/test-data/fastq/1508-1_S1_L004_R2_001.fa
```
</details>


All output is in a `results` folder.

<details><summary>See explaination of <b>results</b> folder</summary>
  
  ```
  results/
    |_ report.html       # detailed breakdown of which processes where run on what input
    |_ timeline.html     # gantt chart-like timeline of each process and how long it ran
    |
    |_ fastqc/           # Contains the html files generated by fastqc quality check
    |_ 0_index/          # Contains the genome index files generated by gatk0
    |_ ....
  ```
  
</details>

-->




