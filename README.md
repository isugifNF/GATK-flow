# GATK

A [Nextflow](https://www.nextflow.io/) wrapper for the [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us), modified from the pipeline described in the Bioinformatic Workbook: [GATK Best Practices Workflow for DNA-Seq](https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gatk-dnaseq-best-practices-workflow.html#gsc.tab=0).

This Nextflow wrapper extends GATK's functionality from DNAseq data, to enabling variant calling from RNA-Seq and long-read sequencing data. It provides a comprehensive solution for diverse sequencing platforms, enhancing variant analysis capabilities.

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

# === Atlas (will need a local install of nextflow and will need the --account "projectname" flag)
module load singularity
NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
```

</details>

```
nextflow run isugifNF/GATK --help
```

<details><summary>See help statement</summary>

```
N E X T F L O W  ~  version 20.10.0
Launching `isugifNF/GATK` [big_kare] - revision: ca139b5b5f
Usage:
   The typical command for running the pipeline is as follows:
   nextflow run main.nf --genome GENOME.fasta --reads "*_{R1,R2}.fastq.gz" -profile slurm,singularity
   nextflow run main.nf --genome GENOME.fasta --reads_file READ_PATHS.txt -profile slurm,singularity

   Mandatory arguments:
    --genome                Genome fasta file, against which reads will be mapped to find SNPs
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
    or
    --genome                Genome fasta file, against which reads will be mapped to find SNPs
    --reads_file            Text file (tab delimited) with three columns [readname left_fastq.gz right_fastq.gz]. Will need full path for files.

    --invariant             Output invariant sites [default:false]
    --long_reads            Long read file in fastq.gz format, will need to specify glob (e.g. "*.fastq.gz")

   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, slurm, singularity, docker [default:local]
    --singularity_img       Singularity image if [-profile singularity] is set [default:'shub://aseetharam/gatk:latest']
    --docker_img            Docker image if [-profile docker] is set [default:'j23414/gatk4']
    --gatk_app              Link to gatk executable [default: 'gatk']
    --bwamem2_app           Link to bwamem2 executable [default: 'bwa-mem2']
    --samtools_app          Link to samtools executable [default: 'samtools']
    --bedtools_app          Link to bedtools executable [default: 'bedtools']
    --datamash_app          Link to datamash executable [default: 'datamash']
    --vcftools_app          Link to vcftools executable [default: 'vcftools']
    
    --star_index_param      Parameters to pass to STAR index [default:'null']
    --star_index_file       Optional: speedup by providing a prebuilt STAR indexed genome [default: 'false']

   Optional other arguments:
    --threads               Threads per process [default:4 for local, 16 for slurm]
    --window                Window size passed to bedtools for gatk [default:100000]
    --queueSize             Maximum jobs to submit to slurm [default:20]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help

NOTE: Genome name should not contain any periods before the ".fasta"
```

</details>

### Singularity Container

Tools required for the workflow are included in the container [aseetharam/gatk:latest](https://github.com/aseetharam/gatk) and should be automatically pulled by nextflow. (Will only need to run `singularity pull` if website connection is unstable.)

<details><summary>If website connection is unstable, pull singularity and use the `-with-singularity` flag</summary>

#### To pull the image

```
singularity pull --name gatk.sif shub://aseetharam/gatk:latest
```

#### Link image to Nextflow using the `-with-singularity` flag.

```
nextflow run isugifNF/GATK \
  --genome "test-data/ref/b73_chr1_150000001-151000000.fasta" \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm \
  -with-singularity gatk.sif
```

</details>

## Test Dataset

A simple DNAseq test dataset (`test-data`) is available on [ISU Box](https://iastate.app.box.com/v/gatk-test-data). This dataset contains a small genome (portion of chr1, B73v5 ), and Illumina short reads for 26 NAM lines (including B73) and B73Ab10 line (27 lines total).
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

This pipeline can process DNAseq, RNAseq and PacBio long read data, which require different arguments and different input files. Please jump to the section based off of your input below.

### GATK variant calling for DNAseq data

Because this pipeline was initially developed to process DNAseq data, we have provided a test DNAseq dataset. You can fetch the pipeline and the test-data folder.

```
# Fetch pipeline
nextflow run isugifNF/GATK --help

# Fetch the test-data folder from ISU box
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
```

The general format of a run with the pipeline is to provide a genome file (`--genome`) and Illumina Paired-End Reads files (`--reads` or `--reads_file`).

```
nextflow run isugifNF/GATK \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm,singularity \
  -resume
```

or

```
nextflow run isugifNF/GATK \
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

#### DNAseq final output

The final output will be in a `results` folder. SNPs will be in the VCF file, probably the file with the longest name (e.g. `first-round_merged_snps-only_snp-only.pass-only.vcf`).

```
results/
  |_ 01_MarkAdapters/            #<= folders contain intermediate files
  |_ 02_MapReads/
  |_ 03_PrepGATK/
  |_ 04_GATK/
  |_ 05_FilterSNPs/
  |  |_ first-round_merged_snps-only_sorted_snp-only.pass-only.vcf     #<= final SNP file
  |
  |_ report.html
  |_ timeline.html               # <= runtime information for all processes

```
### GATK variant calling for RNAseq

In order to process RNAseq, we use the `STAR` aligner. To use the RNAseq variant calling pipeline, make the following parameter changes. 

```
nextflow run isugifNF/GATK \
  --seq "rna" \
  --genome "genome/genome.fa" \
  --reads "data/*.{r1,r2}.fq.gz" \
  --gtf "genome/genome.gtf" \
  -profile slurm,singularity \
  -resume
```

#### RNAseq final output

The final output will be in a `results` folder. SNPs will be in the VCF file, probably the file with the longest name (e.g. `first-round_merged_snps-only_snp-only.pass-only.vcf`).

```
results/
  |_ 01_MarkAdapters/            #<= folders contain intermediate files
  |_ 02_MapReads/
  |_ 03_PrepGATK/
  |_ 04_GATK/
  |_ 05_FilterSNPs/
  |  |_ first-round_merged_snps-only_sorted_snp-only.pass-only.vcf     #<= final SNP file
  |
  |_ report.html
  |_ timeline.html               # <= runtime information for all processes

```

### GATK variant calling for long read data

In order to process long read data obtained from the PacBio platform, we use the `pbmm2` aligner. To use the long read variant calling pipeline, make the following parameter changes. 

```
nextflow run isugifNF/GATK \
  --seq "longread" \
  --genome "genome/genome.fa" \
  --reads "data/long-reads.fastq" \
  -profile slurm,singularity \
  -resume
```

### Long read final output

The final output will be in a `results` folder. SNPs will be in the VCF file, probably the file with the longest name (e.g. `output_snp-only.pass-only.vcf`).

```
results/
  |_ 02_MapReads/
  |_ 03_PrepGATK/
  |_ 04_GATK/
  |_ 05_FilterSNPs/
  |  |_ output_snp-only.pass-only.vcf     #<= final SNP file
  |
  |_ report.html
  |_ timeline.html               # <= runtime information for all processes

```

## Example Runs

Some example runs provided to show nextflow output.

### Example runs for DNAseq data

<details><summary>See example run on <b>Ceres HPC</b> - last update: 14 April 2021</summary>

Runtime: 1 hour 9 minutes and 27 seconds.

```
$ module load nextflow
$ nextflow run main.nf \
  --genome "test-data/ref/b73_chr1_150000001-151000000.fasta" \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  --queueSize 25 \
  -profile slurm,singularity \
  -resume

N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [exotic_poincare] - revision: ca139b5b5f
executor >  slurm (155)
[8c/e6342a] process > FastqToSam (BioSample26)       [100%] 27 of 27 ✔
[6d/13aa51] process > MarkIlluminaAdapters (27_Bi... [100%] 27 of 27 ✔
[dc/032273] process > SamToFastq (20_BioSample24_... [100%] 27 of 27 ✔
[0e/a8d488] process > bwamem2_index (b73_chr1_150... [100%] 1 of 1 ✔
[2f/b3746a] process > bwamem2_mem (20_BioSample24)   [100%] 27 of 27 ✔
[bc/8e43cc] process > CreateSequenceDictionary (b... [100%] 1 of 1 ✔
[fb/c2b500] process > samtools_faidx (b73_chr1_15... [100%] 1 of 1 ✔
[64/b1241a] process > MergeBamAlignment (20_BioSa... [100%] 27 of 27 ✔
[16/ea1c05] process > bedtools_makewindows (b73_c... [100%] 1 of 1 ✔
[34/14a2e1] process > gatk_HaplotypeCaller (chr1:... [100%] 10 of 10 ✔
[a2/6d9b6d] process > merge_vcf                      [100%] 1 of 1 ✔
[dd/187a85] process > vcftools_snp_only (first-ro... [100%] 1 of 1 ✔
[cc/367a9a] process > SortVcf (first-round_merged... [100%] 1 of 1 ✔
[ef/102130] process > calc_DPvalue (first-round_m... [100%] 1 of 1 ✔
[8f/281023] process > VariantFiltration (first-ro... [100%] 1 of 1 ✔
[fc/5d8e7c] process > keep_only_pass (first-round... [100%] 1 of 1 ✔
Completed at: 14-Apr-2021 01:51:21
Duration    : 1h 9m 27s
CPU hours   : 7.3
Succeeded   : 155
```

</details>

<details><summary>See example run on <b>Atlas HPC</b> - last update: 14 April 2021</summary>

Runtime: 50 minutes and 50 seconds.

```
$ module load singularity
$ NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
$ $NEXTFLOW run main.nf \
  --genome "test-data/ref/b73_chr1_150000001-151000000.fasta" \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  --queueSize 50 \
  --account isu_gif_vrsc \
  -profile slurm,singularity \
  -resume

N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [tiny_cori] - revision: ca139b5b5f
executor >  slurm (155)
[92/cfaf35] process > FastqToSam (BioSample24)       [100%] 27 of 27 ✔
[59/37e2e8] process > MarkIlluminaAdapters (20_Bi... [100%] 27 of 27 ✔
[b3/cc67d6] process > SamToFastq (20_BioSample24_... [100%] 27 of 27 ✔
[46/9d4a5f] process > bwamem2_index (b73_chr1_150... [100%] 1 of 1 ✔
[e4/ad114c] process > bwamem2_mem (20_BioSample24)   [100%] 27 of 27 ✔
[a7/044560] process > CreateSequenceDictionary (b... [100%] 1 of 1 ✔
[50/e81af8] process > samtools_faidx (b73_chr1_15... [100%] 1 of 1 ✔
[67/1c03a5] process > MergeBamAlignment (20_BioSa... [100%] 27 of 27 ✔
[f6/d94e63] process > bedtools_makewindows (b73_c... [100%] 1 of 1 ✔
[49/3a6d6c] process > gatk_HaplotypeCaller (chr1:... [100%] 10 of 10 ✔
[68/e88beb] process > merge_vcf                      [100%] 1 of 1 ✔
[55/66c4cf] process > vcftools_snp_only (first-ro... [100%] 1 of 1 ✔
[92/7bad2f] process > SortVcf (first-round_merged... [100%] 1 of 1 ✔
[0b/e824cf] process > calc_DPvalue (first-round_m... [100%] 1 of 1 ✔
[a9/33a4f2] process > VariantFiltration (first-ro... [100%] 1 of 1 ✔
[ef/769ed6] process > keep_only_pass (first-round... [100%] 1 of 1 ✔
Completed at: 14-Apr-2021 01:03:12
Duration    : 50m 50s
CPU hours   : 6.4
Succeeded   : 155
```

</details>

<details><summary>See example run on <b>Nova HPC</b> - last update: 14 April 2021</summary>

Runtime: 1 hour 46 minutes and 17 seconds.

```
$ module load gcc/7.3.0-xegsmw4 nextflow
$ module load singularity
$ nextflow run main.nf \
  --genome "test-data/ref/b73_chr1_150000001-151000000.fasta" \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  --queueSize 25 \
  -profile slurm,singularity \
  -resume

N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [condescending_monod] - revision: ca139b5b5f
executor >  slurm (155)
[5c/b08536] process > FastqToSam (BioSample04)       [100%] 27 of 27 ✔
[88/46d321] process > MarkIlluminaAdapters (27_Bi... [100%] 27 of 27 ✔
[96/200ac5] process > SamToFastq (21_BioSample24_... [100%] 27 of 27 ✔
[4c/8735b5] process > bwamem2_index (b73_chr1_150... [100%] 1 of 1 ✔
[86/06740e] process > bwamem2_mem (21_BioSample24)   [100%] 27 of 27 ✔
[c0/3e6521] process > CreateSequenceDictionary (b... [100%] 1 of 1 ✔
[e2/856737] process > samtools_faidx (b73_chr1_15... [100%] 1 of 1 ✔
[25/529408] process > MergeBamAlignment (21_BioSa... [100%] 27 of 27 ✔
[40/ca5ca7] process > bedtools_makewindows (b73_c... [100%] 1 of 1 ✔
[8a/0d6f00] process > gatk_HaplotypeCaller (chr1:... [100%] 10 of 10 ✔
[96/0b957b] process > merge_vcf                      [100%] 1 of 1 ✔
[96/1f9848] process > vcftools_snp_only (first-ro... [100%] 1 of 1 ✔
[60/edb33d] process > SortVcf (first-round_merged... [100%] 1 of 1 ✔
[b0/372a4e] process > calc_DPvalue (first-round_m... [100%] 1 of 1 ✔
[f3/cfb966] process > VariantFiltration (first-ro... [100%] 1 of 1 ✔
[86/024c8b] process > keep_only_pass (first-round... [100%] 1 of 1 ✔
Completed at: 14-Apr-2021 02:17:48
Duration    : 1h 46m 17s
CPU hours   : 6.6
Succeeded   : 155
```

</details>

### Example run for RNAseq data

<details><summary>see example run on <b>NOVA hpc</b> - last update: 8 august 2023</summary>

runtime: 4 hours 38 minutes and 54 seconds.

```
$ module load gcc/7.3.0-xegsmw4 nextflow
$ module load singularity
$ nextflow run main.nf \
  --seq "rna" \
  --genome "genome/s_lycopersicum_chromosomes.3.00.fa" \
  --reads "01_data/s*.{r1,r2}.fq.gz" \
  --queuesize 25 \
  --java_options "-xmx80g -xx:+useparallelgc -djava.io.tmpdir=/work/gif3/satheesh/2023_gatk_ms/tmp" \
  --gtf "genome/itag3.0_gene_models.gtf" \
  -profile slurm \
  -resume

N E X T F L O W  ~  version 22.10.4
launching `/work/gif3/satheesh/programs/gatk_testing/gatk/main.nf` [chaotic_brahmagupta] dsl2 - revision: d07ccde0d9
executor >  slurm (8304)
[c4/48a5ae] process > fastqtosam (subsampled_p7r1)   [100%] 1 of 1 ✔
[ee/c542a3] process > markilluminaadapters (1_sub... [100%] 1 of 1 ✔
[f8/a5cd97] process > samtofastq_rna (1_subsample... [100%] 1 of 1 ✔
[8b/c7dafa] process > star_index (s_lycopersicum_... [100%] 1 of 1 ✔
[2c/f96408] process > star_align (1_subsampled_p7r1) [100%] 1 of 1 ✔
[ee/72b0af] process > createsequencedictionary (s... [100%] 1 of 1 ✔
[29/6454c2] process > samtools_faidx (s_lycopersi... [100%] 1 of 1 ✔
[af/2ff081] process > mergebamalignment_rna (1_su... [100%] 1 of 1 ✔
[51/589571] process > markduplicates (1_subsample... [100%] 1 of 1 ✔
[4c/53bb96] process > splitncigarreads (1_subsamp... [100%] 1 of 1 ✔
[ee/79c9ea] process > bedtools_makewindows (s_lyc... [100%] 1 of 1 ✔
[f6/a83bb1] process > gatk_haplotypecaller_rna (s... [100%] 8287 of 8287, fai...
[5c/efa2eb] process > merge_vcf                      [100%] 1 of 1 ✔
[ec/6b40be] process > vcftools_snp_only (first_ro... [100%] 1 of 1 ✔
[ab/604154] process > sortvcf (first_round_merged... [100%] 1 of 1 ✔
[28/765be4] process > calc_dpvalue (first_round_m... [100%] 1 of 1 ✔
[5b/0b530f] process > variantfiltration (first_ro... [100%] 1 of 1 ✔
[99/878369] process > keep_only_pass (first_round... [100%] 1 of 1 ✔
warn: failed to render execution report -- see the log file for details
warn: failed to render execution timeline -- see the log file for details
completed at: 08-aug-2023 21:37:02
duration    : 4h 38m 54s
cpu hours   : 6.3
succeeded   : 8'303
```



</details>

### Example run for long read data

<details><summary>see example run on <b>NOVA hpc</b> - last update: 31 October 2023</summary>

runtime: 41 minutes and 23 seconds.

```
$ module load gcc/7.3.0-xegsmw4 nextflow
$ module load singularity
$ nextflow run main.nf \
  --seq "longread" \
  --genome "genome/ecoli_2000001_3000000.fasta" \
  --reads "01_Data/CP007799_2000001-3000000.fastq" \
  --queueSize 25 \
  --java_options "-Xmx80g -XX:+UseParallelGC -Djava.io.tmpdir=/work/gif3/satheesh/2023_GATK_ms/tmp" \
  -profile slurm,singularity \
  -resume

N E X T F L O W  ~  version 22.10.4
launching `/work/gif3/satheesh/programs/gatk_testing/gatk/main.nf` [chaotic_brahmagupta] dsl2 - revision: d07ccde0d9
executor >  slurm (21)
[c8/aa3601] process > LONGREAD_VARIANT_CALLING:Cr... [100%] 1 of 1 ✔
[51/b76b1a] process > LONGREAD_VARIANT_CALLING:sa... [100%] 1 of 1 ✔
[be/0cb73f] process > LONGREAD_VARIANT_CALLING:be... [100%] 1 of 1 ✔
[26/4f4559] process > LONGREAD_VARIANT_CALLING:pb... [100%] 1 of 1 ✔
[d1/78ce63] process > LONGREAD_VARIANT_CALLING:pb... [100%] 1 of 1 ✔
[bb/1e7433] process > LONGREAD_VARIANT_CALLING:Ma... [100%] 1 of 1 ✔
[aa/4c9f99] process > LONGREAD_VARIANT_CALLING:ga... [100%] 10 of 10 ✔
[a6/142af0] process > LONGREAD_VARIANT_CALLING:Co... [100%] 1 of 1 ✔
[fa/069ebb] process > LONGREAD_VARIANT_CALLING:Ge... [100%] 1 of 1 ✔
[09/ea1b74] process > LONGREAD_VARIANT_CALLING:ca... [100%] 1 of 1 ✔
[65/cc76dd] process > LONGREAD_VARIANT_CALLING:Va... [100%] 1 of 1 ✔
[77/452275] process > LONGREAD_VARIANT_CALLING:ke... [100%] 1 of 1 ✔
Completed at: 31-Oct-2023 14:39:44
Duration    : 41m 23s
CPU hours   : 5.7
Succeeded   : 21
```

</details>


