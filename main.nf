#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { DNA_VARIANT_CALLING } from './subworkflows/local/dna_variant_calling/main.nf'

include { RNA_VARIANT_CALLING } from './subworkflows/local/rna_variant_calling/main.nf'

def helpMsg() {
  log.info """
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

    --invariant             Output invariant sites [default:false]

   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, slurm, singularity, docker [default:local]
    --singularity_img       Singularity image if [-profile singularity] is set [default:'${params.singularity_img}']
    --docker_img            Docker image if [-profile docker] is set [default:'${params.docker_img}']
    --gatk_app              Link to gatk executable [default: '$gatk_app']
    --bwamem2_app           Link to bwamem2 executable [default: '$bwamem2_app']
    --samtools_app          Link to samtools executable [default: '$samtools_app']
    --bedtools_app          Link to bedtools executable [default: '$bedtools_app']
    --datamash_app          Link to datamash executable [default: '$datamash_app']
    --vcftools_app          Link to vcftools executable [default: '$vcftools_app']

    --star_index_param      Parameters to pass to STAR index [default:'${params.star_index_param}']
    --star_index_file       Optional: speedup by providing a prebuilt STAR indexed genome [default: '${params.star_index_file}']

   Optional other arguments:
    --java_options          Java options for gatk [default:'${java_options}']
    --threads               Threads per process [default:4 for local, 16 for slurm]
    --window                Window size passed to bedtools for gatk [default:${params.window}]
    --queueSize             Maximum jobs to submit to slurm [default:${params.queueSize}}]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help

"""
}

if(params.help){
  helpMsg()
  exit 0
}

def parameters_valid = ['help','outdir',
  'genome','gtf','reads','reads_file','invariant','seq',
  'singularity_img','docker_img',
  'gatk_app','star_app','star_index_params','star_index_file','bwamem2_app','samtools_app','bedtools_app','datamash_app','vcftools_app',
  'java_options','window','queueSize','queue-size','account', 'threads'] as Set

def parameter_diff = params.keySet() - parameters_valid
if (parameter_diff.size() != 0){
   exit 1, "[Pipeline error] Parameter(s) $parameter_diff is(are) not valid in the pipeline!\n"
}

workflow {
  // == Read in genome and reads channels
  if(params.genome) {
    genome_ch = channel.fromPath(params.genome, checkIfExists:true)
      | view {file -> "Genome file : $file "}
  } else {
    exit 1, "[Missing File(s) Error] This pipeline requires a reference '--genome [GENOME.fasta]' \n"
  }

  if (params.reads) {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
      | view {files -> "Read files : $files "}
  } else if (params.reads_file) {
    reads_ch = channel.fromPath(params.reads_file, checkIfExists:true)
      | splitCsv(sep:'\t')
      | map { n -> [ n.getAt(0), [n.getAt(1), n.getAt(2)]] }
      | view {files -> "Read files : $files "}
  } else {
    exit 1, "[Missing File(s) Error] This pipeline requires either paired-end read files as a glob '--reads [*_{r1,r2}.fq.gz]' or as a tab-delimited text file '--reads_file [READS_FILE.txt]'\n"
  }

  if (params.seq == "dna"){
    DNA_VARIANT_CALLING(genome_ch, reads_ch)
  } else if (params.seq == "rna") {
    gtf_ch = channel.fromPath(params.gtf, checkIfExists:true)
    RNA_VARIANT_CALLING(genome_ch, reads_ch, gtf_ch)    
  }
}
