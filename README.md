# Maize WGS Build

<hr/>

A [Nextflow](https://www.nextflow.io/) wrapper for the [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us), modified from the pipeline described in the Bioinformatic Workbook: [GATK Best Practices Workflow for DNA-Seq](https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gatk-dnaseq-best-practices-workflow.html#gsc.tab=0).

<!-- The benefits of Nextflow include:

* write once, run anywhere (`configs/*.config` for singularity, slurm, local)
* checkpointing and runtime reports
* customizing for a particular HPC
-->

<!--### Dependencies

For portability, the dependencies for the GATK pipeline are provided as a singularity image. To avoid singularity, the individual programs (`bwa`, `samtools`, `picard`, `bedtools`, `gatk`, `vcftools`) can be configured directly, see help statement (e.g. `--samtools_app`) under "Installation".

* [Nextflow](https://www.nextflow.io/)
* Singularity Image File
* Input files:
  * genome file (`some_genome.fasta`) 
  * Illumina Paired End reads (`reads_R1.fastq`, `reads_R2.fastq.gz`)
-->

### Installation

You will need a working version of nextflow, [see here](https://www.nextflow.io/docs/latest/getstarted.html#requirements) on how to install nextflow. Nextflow modules are avialble on some of the HPC computing resources.

<details><summary>See modules on HPC clusters</summary>

```
# === Nova
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Condo
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Ceres
module load nextflow
module load singularity
NEXTFLOW=nextflow

# === Atlas
module load singularity
NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
```

</details>


```
git clone https://github.com/HuffordLab/Maize_WGS_Build.git
cd Maize_WGS_Build

nextflow run main.nf --help
```

<details><summary>See help statement</summary>

```
N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [zen_woese] - revision: 0516af2de3
Usage:
   The typical command for running the pipeline is as follows:
   nextflow run main.nf --genome GENOME.fasta --reads "*_{R1,R2}.fastq.gz" -profile singularity
   nextflow run main.nf --genome GENOME.fasta --reads_file READ_PATHS.txt -profile singularity

   Mandatory arguments:
    --genome                Genome fasta file, against which reads will be mapped to find SNPs
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
    or
    --genome                Genome fasta file, against which reads will be mapped to find SNPs
    --reads_file            Text file (tab delimited) with three columns [readname left_fastq.gz right_fastq.gz]. Will need full path for files.

   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, condo, atlas, singularity [default:local]
    --singularity_img       Singularity image if [-profile singularity] is set [default:'shub://aseetharam/gatk:latest']
    --bwa_app               Link to bwa executable [default: 'bwa']
    --samtools_app          Link to samtools executable [default: 'samtools']
    --picard_app            Link to picard executable [default: 'picard'], might want to change to "java -jar ~/PICARD_HOME/picard.jar"
    --bedtools_app          Link to bedtools executable [default: 'bedtools']
    --gatk_app              Link to gatk executable [default: 'gatk']
    --datamash_app          Link to datamash executable [default: 'datamash']
    --vcftools_app          Link to vcftools executable [default: 'vcftools']

   Optional other arguments:
    --window                Window size passed to bedtools for gatk [default:100000]
    --queueSize             Maximum jobs to submit to slurm [default:18]
    --help
```

</details>

### Test Dataset

A simple test dataset (`test-data`) is available on [ISU Box](https://iastate.app.box.com/v/gatk-test-data). This dataset contains a small genome (portion of chr1, B73v5 ), and Illumina short reads for 26 NAM lines (including B73) and B73Ab10 line (27 lines total).
Only the reads that map to the region of the v5 genome is included, so that this can be tested quickly.
There are examples of multiple files belonging to same NAM line as well as single file per NAM line to make sure both conditions works correctly.
The end VCF file should have exactly 27 individuals (lines) in them.


```
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
```

<!-- ### Container

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
--> 

### Running the Pipeline

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

Fetch the pipeline and fetch the test-data folder.

```
# Fetch repo
git clone https://github.com/HuffordLab/Maize_WGS_Build.git
cd Maize_WGS_Build

# Fetch the test-data folder from ISU box
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
```

If you are on a HPC (Nova/Condo/Ceres/Atlas), it is highly recommend to use the `submit_nf.slurm` script. The `# === Load Modules` section will need to be modified to get nextflow and singulariy running.

<details><summary>See Module Changes</summary>

```
# === Nova
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Condo
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Ceres
module load nextflow
module load singularity
NEXTFLOW=nextflow

# === Atlas
module load singularity
NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
```

</details>

  
 Example run on Atlas with 27 Illumina paired-end reads (listed in `my_group.txt`) against genome (`ref/b73_chr1_150000001-151000000.fasta`).
 
  ```
  nextflow run HuffordLab/Maize_WGS_Build \
    -profile atlas,singularity \
    --reads_file my_group.txt \
    --genome test-data/ref/b73_chr1_150000001-151000000.fasta
    
  executor >  slurm (156)
  [b9/51a78c] process > prep_genome:fasta_sort (b73... [100%] 1 of 1 ✔
  [ae/b743bd] process > prep_genome:fasta_bwa_index... [100%] 1 of 1 ✔
  [f5/4e914b] process > prep_genome:fasta_samtools_... [100%] 1 of 1 ✔
  [be/d21a0d] process > prep_genome:fasta_picard_di... [100%] 1 of 1 ✔
  [fc/b54240] process > prep_reads:paired_FastqToSA... [100%] 27 of 27 ✔
  [ad/50c6b6] process > prep_reads:BAM_MarkIllumina... [100%] 27 of 27 ✔
  [48/0721fa] process > map_reads:BAM_SamToFastq (B... [100%] 27 of 27 ✔
  [82/7aa82d] process > map_reads:run_bwa_mem (B24_... [100%] 27 of 27 ✔
  [96/56f2c6] process > run_MergeBamAlignment (B02)    [100%] 27 of 27 ✔
  [82/c5cc00] process > fai_bedtools_makewindows (b... [100%] 1 of 1 ✔
  [a4/9badc8] process > run_gatk_snp (chr1:900001-9... [100%] 10 of 10 ✔
  [49/17f481] process > merge_vcf                      [100%] 1 of 1 ✔
  [02/c5fe0a] process > vcftools_snp_only (first-ro... [100%] 1 of 1 ✔
  [90/312c78] process > run_SortVCF (first-round_me... [100%] 1 of 1 ✔
  [5c/fd7ee1] process > calc_DPvalue (first-round_m... [100%] 1 of 1 ✔
  [10/26a37d] process > gatk_VariantFiltration (fir... [100%] 1 of 1 ✔
  [d5/d4fad2] process > keep_only_pass (first-round... [100%] 1 of 1 ✔
  2260.74

  /project/isu_gif_vrsc/Jennifer/github/blank/work/d5/d4fad26c33b234ab856e600353ebb0/first-round_merged_snps-only_snp-only.pass-only.vcf
  Completed at: 15-Oct-2020 14:00:31
  Duration    : 51m 8s
  CPU hours   : 6.1
  Succeeded   : 156
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
