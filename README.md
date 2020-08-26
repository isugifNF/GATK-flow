# Maize_WGS_Build
Code and Data to include in WGS Build

Start with the workflow documented on [Bioinformatic Workbook GATK DNAseq Best Practices](https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gatk-dnaseq-best-practices-workflow.html#gsc.tab=0)

### Running the pipeline

If on a local laptop with nextflow installed:

```
nextflow run HuffordLab/Maize_WGS_Build
```

If on HPCC Condo:

```
module load gcc/7.3.0-xegsmw4 nextflow
nextflow run HuffordLab/Maize_WGS_Build -profile condo
```


### Test dataset

A simple test dataset is available (here)[/test-data]. This dataset contains a small genome (portion of chr1, B73v5), and Illumina short reads for 26 NAM lines (including B73) and B73Ab10 line (27 lines total).
Only the reads that map to the region of the v5 genome is included, so that this can be tested quickly.
There are examples of multiple files belonging to same NAM line as well as single file per NAM line to make sure both conditions works correctly.
The end VCF file should have exactly 27 individuals (lines) in them.

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
singularity exec gatk.sif java -jar $GATKHOME/$GATK
singularity exec gatk.sif java -jar $PICARDHOME/picard.jar
```

