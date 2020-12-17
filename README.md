# GATK


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

## Installation

You will need a working version of nextflow, [see here](https://www.nextflow.io/docs/latest/getstarted.html#requirements) on how to install nextflow. Nextflow modules are avialable on some of the HPC computing resources.

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
# singularity already available, no need for module
NEXTFLOW=nextflow

# === Atlas
module load singularity
NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
```

</details>

<!--
```
git clone https://github.com/HuffordLab/Maize_WGS_Build.git
cd Maize_WGS_Build

nextflow run main.nf --help
```
-->

```
git clone https://github.com/isugfNF/GATK.git
cd GATK

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
    --account               HPC account name for slurm sbatch, atlas and ceres may require this
    --help
```

</details>

### Singularity Container

Tools required for the workflow are included in the container [aseetharam/gatk:latest](https://github.com/aseetharam/gatk) and should be automatically pulled by nextflow. (Will only need to run `singularity pull` if website connection is unstable.)

<details><summary>More Info</summary>

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

</details>


## Test Dataset

A simple test dataset (`test-data`) is available on [ISU Box](https://iastate.app.box.com/v/gatk-test-data). This dataset contains a small genome (portion of chr1, B73v5 ), and Illumina short reads for 26 NAM lines (including B73) and B73Ab10 line (27 lines total).
Only the reads that map to the region of the v5 genome is included, so that this can be tested quickly.
There are examples of multiple files belonging to same NAM line as well as single file per NAM line to make sure both conditions works correctly.
The end VCF file should have exactly 27 individuals (lines) in them.


```
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz

ls -1 test-data/
#> fastq          # <= folder of paired reads
#> read-group.txt
#> ref            # <= folder containing one genome reference
```



## Running the Pipeline

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

<!--
```
# Fetch repo
git clone https://github.com/HuffordLab/Maize_WGS_Build.git
cd Maize_WGS_Build

# Fetch the test-data folder from ISU box
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
```
-->

```
# Fetch repo
git clone https://github.com/isugifNF/GATK.git
cd GATK

# Fetch the test-data folder from ISU box
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
```

The general format of a run with the pipeline is to provide a genome file (`--genome`) and Illumina Paired-End Reads files (`--reads` or `--reads_file`).

```
nextflow run main.nf \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm,singularity \
  -resume
```

or 

```
nextflow run main.nf \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads_file read-path.txt \
  -profile slurm,singularity \
  -resume
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
# singularity already available, no need for module
NEXTFLOW=nextflow

# === Atlas
module load singularity
NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
```

</details>


### Using the --reads_file

Instead of using a pattern to specify paired reads (`"test-data/fastq/*_{R1,R2}.fastq.gz"`), we can use a tab-delimited textfile to specify the path to left and right files. The textfile should contain three columns: readname, left read path, right read path.

In our case, for the `test-data` run the following one-liner to generate a `read-path.txt`.

```
for f in test-data/fastq/*_R1.fastq.gz; do echo -e "$(basename $f |cut -f 1 -d "_")\t$(realpath $f)\t$(realpath $f | sed 's/_R1.fastq.gz/_R2.fastq.gz/g')"; done > read-path.txt
```

<details><summary>See example <b>read-path.txt</b></summary>

```
BioSample01	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample01_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample01_R2.fastq.gz
BioSample02	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample02_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample02_R2.fastq.gz
BioSample03	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample03_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample03_R2.fastq.gz
BioSample04	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample04_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample04_R2.fastq.gz
BioSample05	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample05_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample05_R2.fastq.gz
```

</details>

### Final Output

The Final output will be in a `results` folder. SNPs will be in the VCF file, probably the file with the longest name (e.g. `first-round_merged_snps-only_snp-only.pass-only.vcf`).

<details><summary>See <b>results</b> folder</summary>

```
ls -ltrh results/
total 5.8M
drwxr-s--- 2 user proj 4.0K Oct 12 23:55 sort_fasta/    # <= folders contain intermediate files
drwxr-s--- 2 user proj 4.0K Oct 12 23:55 samtools/
drwxr-s--- 2 user proj 4.0K Oct 12 23:55 bwa/
drwxr-s--- 2 user proj 4.0K Oct 12 23:56 bedtools/
drwxr-s--- 2 user proj 4.0K Oct 12 23:56 bwa_mem/
drwxr-s--- 2 user proj 4.0K Oct 13 00:01 gatk/
drwxr-s--- 2 user proj 4.0K Oct 13 00:01 vcftools/
drwxr-s--- 2 user proj 4.0K Oct 13 00:01 picard/
lrwxrwxrwx 1 user proj  132 Oct 13 00:02 first-round_merged_snps-only.marked.vcf
lrwxrwxrwx 1 user proj  144 Oct 13 00:02 first-round_merged_snps-only_snp-only.pass-only.vcf # <= Final SNP file
-rw-r----- 1 user proj  16K Oct 13 00:02 timeline.html  # <= shows runtime for each portion
-rw-r----- 1 user proj 2.9M Oct 13 00:02 report.html    # <= shows resource use
```

</details>

## Example Runs

Some example runs provided to show nextflow output. These were run using the earlier forked version `HuffordLab/Maize_WGS_Build`

<details><summary>See example run on <b>Ceres HPC</b></summary>

Runtime: 1 hour 7 minutes and 17 seconds.

```
$ nextflow run main.nf \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm,singularity \
  -resume

N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [extravagant_sinoussi] - revision: d5f8cdb041
WARN: It appears you have never run this project before -- Option `-resume` is ignored
executor >  slurm (156)
[58/16126c] process > prep_genome:fasta_sort (b73... [100%] 1 of 1 ✔
[a8/a73011] process > prep_genome:fasta_bwa_index... [100%] 1 of 1 ✔
[43/406cf5] process > prep_genome:fasta_samtools_... [100%] 1 of 1 ✔
[27/647551] process > prep_genome:fasta_picard_di... [100%] 1 of 1 ✔
[05/dc6a6a] process > prep_reads:paired_FastqToSA... [100%] 27 of 27 ✔
[18/c53c3d] process > prep_reads:BAM_MarkIllumina... [100%] 27 of 27 ✔
[9d/f099dd] process > map_reads:BAM_SamToFastq (B... [100%] 27 of 27 ✔
[5a/95f9dd] process > map_reads:run_bwa_mem (BioS... [100%] 27 of 27 ✔
[6a/56f359] process > run_MergeBamAlignment (BioS... [100%] 27 of 27 ✔
[24/9d4e1d] process > fai_bedtools_makewindows (b... [100%] 1 of 1 ✔
[f5/3a7b48] process > run_gatk_snp (chr1:900001-9... [100%] 10 of 10 ✔
[f8/9dd8e8] process > merge_vcf                      [100%] 1 of 1 ✔
[40/0b74b8] process > vcftools_snp_only (first-ro... [100%] 1 of 1 ✔
[a6/9c6152] process > run_SortVCF (first-round_me... [100%] 1 of 1 ✔
[af/2f9784] process > calc_DPvalue (first-round_m... [100%] 1 of 1 ✔
[bf/c7d576] process > gatk_VariantFiltration (fir... [100%] 1 of 1 ✔
[51/d58598] process > keep_only_pass (first-round... [100%] 1 of 1 ✔
2265.1

/lustre/project/isu_gif_vrsc/jenchang/_wrkspc/2020-11-15_NF_account/Maize_WGS_Build/work/51/d585988b0193a9cf0aceb653e468de/first-round_merged_snps-only_snp-only.pass-only.vcf
Completed at: 15-Nov-2020 15:24:53
Duration    : 1h 7m 17s
CPU hours   : 7.6
Succeeded   : 156

```

</details>

<details><summary>See example run on <b>Atlas HPC</b></summary>

Example run on Atlas with 27 Illumina paired-end reads (listed in `my_group.txt`) against genome (`ref/b73_chr1_150000001-151000000.fasta`).

Runtime: 50 minutes and 7 seconds.
 
```
$ nextflow run main.nf \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile atlas,singularity \
  -resume
  
N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [awesome_poincare] - revision: d5f8cdb041 
executor >  slurm (156)
[b6/77467c] process > prep_genome:fasta_sort (b73... [100%] 1 of 1 ✔
[39/4f61d2] process > prep_genome:fasta_bwa_index... [100%] 1 of 1 ✔
[40/c4d2df] process > prep_genome:fasta_samtools_... [100%] 1 of 1 ✔
[d3/2c6941] process > prep_genome:fasta_picard_di... [100%] 1 of 1 ✔
[f4/25d924] process > prep_reads:paired_FastqToSA... [100%] 27 of 27 ✔
[b2/693604] process > prep_reads:BAM_MarkIllumina... [100%] 27 of 27 ✔
[a8/242f07] process > map_reads:BAM_SamToFastq (B... [100%] 27 of 27 ✔
[95/afae47] process > map_reads:run_bwa_mem (BioS... [100%] 27 of 27 ✔
[f6/4425e9] process > run_MergeBamAlignment (BioS... [100%] 27 of 27 ✔
[73/eebdca] process > fai_bedtools_makewindows (b... [100%] 1 of 1 ✔
[fe/b96a27] process > run_gatk_snp (chr1:900001-9... [100%] 10 of 10 ✔
[8b/64c2b5] process > merge_vcf                      [100%] 1 of 1 ✔
[48/ab3ec5] process > vcftools_snp_only (first-ro... [100%] 1 of 1 ✔
[0d/2851e0] process > run_SortVCF (first-round_me... [100%] 1 of 1 ✔
[d9/6eff09] process > calc_DPvalue (first-round_m... [100%] 1 of 1 ✔
[ec/254329] process > gatk_VariantFiltration (fir... [100%] 1 of 1 ✔
[20/b0b7ae] process > keep_only_pass (first-round... [100%] 1 of 1 ✔
2265.1

/project/isu_gif_vrsc/Jennifer/github/Maize_WGS_Build/work/20/b0b7ae92102279447059de61707faa/first-round_merged_snps-only_snp-only.pass-only.vcf
Completed at: 28-Oct-2020 13:34:41
Duration    : 50m 7s
CPU hours   : 6.2
Succeeded   : 156   
```
  
</details>

<details><summary>See example run on <b>Condo HPC</b></summary>

Runtime: 2 hours 5 minutes and 39 seconds.

```
$ nextflow run main.nf \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm,singularity \
  -resume
  
N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [clever_monod] - revision: d5f8cdb041
WARN: It appears you have never run this project before -- Option `-resume` is ignored
executor >  slurm (156)
[69/9f3959] process > prep_genome:fasta_sort (b73... [100%] 1 of 1 ✔
[f9/c3d116] process > prep_genome:fasta_bwa_index... [100%] 1 of 1 ✔
[b8/d79a0e] process > prep_genome:fasta_samtools_... [100%] 1 of 1 ✔
[22/9ebdcb] process > prep_genome:fasta_picard_di... [100%] 1 of 1 ✔
[a3/3b449a] process > prep_reads:paired_FastqToSA... [100%] 27 of 27 ✔
[6b/52d8e1] process > prep_reads:BAM_MarkIllumina... [100%] 27 of 27 ✔
[d4/bebcc3] process > map_reads:BAM_SamToFastq (B... [100%] 27 of 27 ✔
[42/367b85] process > map_reads:run_bwa_mem (BioS... [100%] 27 of 27 ✔
[ca/71fa06] process > run_MergeBamAlignment (BioS... [100%] 27 of 27 ✔
[ea/a3adc0] process > fai_bedtools_makewindows (b... [100%] 1 of 1 ✔
[f4/683387] process > run_gatk_snp (chr1:900001-9... [100%] 10 of 10 ✔
[45/bc1e85] process > merge_vcf                      [100%] 1 of 1 ✔
[f4/5e9035] process > vcftools_snp_only (first-ro... [100%] 1 of 1 ✔
[2d/58f2c9] process > run_SortVCF (first-round_me... [100%] 1 of 1 ✔
[df/c75b2a] process > calc_DPvalue (first-round_m... [100%] 1 of 1 ✔
[3c/9cec07] process > gatk_VariantFiltration (fir... [100%] 1 of 1 ✔
[90/6b176a] process > keep_only_pass (first-round... [100%] 1 of 1 ✔
/work/GIF/jenchang/github/Maize_WGS_Build/work/90/6b176a1ae430c87dca4745359652e3/first-round_merged_snps-only_snp-only.pass-only.vcf
Completed at: 28-Oct-2020 15:03:08
Duration    : 2h 5m 39s
CPU hours   : 10.1
Succeeded   : 156
```

</details>

<details><summary>See example run on <b>MacOS</b> laptop where dependencies are locally installed</summary>

(2) On MacOS laptop where dependencies are locally installed:

Need to rerun the local version... will be updated.

```
$ nextflow run main.nf \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  --picard_app "java -jar ~/bin/picard.jar" \
  -profile local

N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [amazing_rubens] - revision: 66f7e69455
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
