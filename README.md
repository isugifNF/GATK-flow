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

(1) On HPCC (condo/nova/atlas)

```
# Fetch repo
git clone https://github.com/HuffordLab/Maize_WGS_Build.git
cd Maize_WGS_Build

# Fetch dataset from ISU box here and place in `test-data` folder
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
```

(1a) Condo HPC, using singularity

  ```
  module load gcc/7.3.0-xegsmw4 nextflow
  module load singularity
  nextflow run main_temp.nf -profile slurm,singularity -resume
  ```
  
(1b) Atlas HPC, using singularity

  ```
  module load singularity
  NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
  ${NEXTFOW} run main_temp.nf -profile atlas,singularity -resume
  ```

(2) On MacOS laptop where dependencies are locally installed:

```
$ nextflow run main_temp.nf --picard_app "java -jar ~/bin/picard.jar"

N E X T F L O W  ~  version 20.07.1
Launching `main_temp.nf` [amazing_rubens] - revision: 66f7e69455
[a4/41e1ad] process > prep_genome:fasta_sort (b73_chr1_150000001-151000000.fasta)                  [100%] 1 of 1, cached: 1 ✔
[f4/b63b1b] process > prep_genome:fasta_bwa_index (b73_chr1_150000001-151000000_sorted.fasta)      [100%] 1 of 1, cached: 1 ✔
[da/436b26] process > prep_genome:fasta_samtools_faidx (b73_chr1_150000001-151000000_sorted.fasta) [100%] 1 of 1, cached: 1 ✔
[18/1d2871] process > prep_genome:fasta_picard_dict (b73_chr1_150000001-151000000_sorted.fasta)    [100%] 1 of 1, cached: 1 ✔
[f8/12a295] process > prep_reads:paired_FastqToSAM (BioSample05)                                   [100%] 3 of 3, cached: 3 ✔
[41/2ca9e5] process > prep_reads:BAM_MarkIlluminaAdapters (BioSample05.bam)                        [100%] 3 of 3, cached: 3 ✔
[35/7b205c] process > map_reads:BAM_SamToFastq (BioSample05_marked.bam)                            [100%] 3 of 3, cached: 3 ✔
[4b/529b7c] process > map_reads:run_bwa_mem (BioSample05_marked_interleaved.fq)                    [100%] 3 of 3, cached: 3 ✔
[d7/4d996a] process > run_MergeBamAlignment (BioSample05)                                          [100%] 3 of 3, cached: 3 ✔
[f4/ebf7ef] process > fai_bedtools_makewindows (b73_chr1_150000001-151000000_sorted.fasta.fai)     [100%] 1 of 1, cached: 1 ✔
[5a/44739d] process > run_gatk_snp (chr1:900001-999999)                                            [100%] 10 of 10, cached: 10 ✔
[5e/ca4c3a] process > merge_vcf                                                                    [100%] 1 of 1, cached: 1 ✔
[4e/43bf37] process > vcftools_snp_only (first-round_merged.vcf)                                   [100%] 1 of 1, cached: 1 ✔
[25/924b29] process > run_SortVCF (1)                                                              [100%] 1 of 1, cached: 1 ✔
[31/b276cb] process > calc_DPvalue (first-round_merged_snps-only.sorted.vcf)                       [100%] 1 of 1, cached: 1 ✔
[4a/17fb0b] process > gatk_VariantFiltration (1)                                                   [100%] 1 of 1, cached: 1 ✔
[88/bcb457] process > keep_only_pass (1)                                                           [100%] 1 of 1, cached: 1 ✔
406.235

/Users/jenchang/Desktop/new/Maize_WGS_Build/work/88/bcb45710b109051ac54bbb0b2fb682/first-round_merged_snps-only_snp-only.pass-only.vcf
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
