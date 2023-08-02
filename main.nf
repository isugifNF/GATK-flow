#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { FastqToSam;
          MarkIlluminaAdapters;
          CreateSequenceDictionary;
          samtools_faidx;
          bedtools_makewindows;
          gatk_HaplotypeCaller;
          gatk_HaplotypeCaller_invariant;
          CombineGVCFs;
          GenotypeGVCFs;
          merge_vcf;
          vcftools_snp_only;
          SortVcf;
          calc_DPvalue;
          VariantFiltration;
          keep_only_pass; } from './modules/GATK.nf'

include { SamToFastq as SamToFastq_RNA;
          STAR_index;
          STAR_align;
          MergeBamAlignment as MergeBamAlignment_RNA; 
          MarkDuplicates; 
          SplitNCigarReads; } from './modules/RNAseq.nf'

include { SamToFastq as SamToFastq_DNA
          bwamem2_index;
          bwamem2_mem; 
          MergeBamAlignment as MergeBamAlignment_DNA; } from './modules/DNAseq.nf'

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
  'gatk_app','star_app','star_index_params','bwamem2_app','samtools_app','bedtools_app','datamash_app','vcftools_app',
  'java_options','window','queueSize','queue-size','account', 'threads'] as Set

def parameter_diff = params.keySet() - parameters_valid
if (parameter_diff.size() != 0){
   exit 1, "[Pipeline error] Parameter(s) $parameter_diff is(are) not valid in the pipeline!\n"
}


// if(!params.genome) {
//   log.info"""
// #===============
//   ERROR: --genome GENOME.fasta    A reference genome file is required!
// #===============
//   """
//   helpMsg()
//   exit 0
// }
// 
// if(!params.reads & !params.reads_file){
//   log.info"""
// #===============
//   ERROR: --reads "*_{r1,r2}.fq.gz"     Paired-end read files are required! Either as a glob or as a tab-delimited text file
//          --reads_file READS_FILE.txt
// #===============
//   """
//   helpMsg()
//   exit 0
// }

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

  // == Since one sample may be run on multiple lanes
  i = 1

  // == Prepare mapped and unmapped read files
  if(params.seq == "dna"){
    cleanreads_ch = reads_ch
    | map { n -> [n.getAt(0), n.getAt(1), "${i++}_"+n.getAt(0)] }
    | FastqToSam
    | MarkIlluminaAdapters
    | SamToFastq_DNA
    | map { n -> [ n.getAt(0).replaceFirst("_marked",""), [ n.getAt(1), n.getAt(2)] ] }
  } else if( params.seq == "rna"){
    cleanreads_ch = reads_ch
    | map { n -> [n.getAt(0), n.getAt(1), "${i++}_"+n.getAt(0)] }
    | FastqToSam
    | MarkIlluminaAdapters
    | SamToFastq_RNA
    | map { n -> [ n.getAt(0).replaceFirst("_marked",""), [ n.getAt(1), n.getAt(2)] ] }
  }

  if(params.seq == "dna"){
    mapped_ch = genome_ch
    | bwamem2_index
    | combine(cleanreads_ch)
    | bwamem2_mem
  } else if( params.seq == "rna"){
    gtf_ch = channel.fromPath(params.gtf, checkIfExists:true)

    mapped_ch = genome_ch
    | combine(gtf_ch)
    | STAR_index
    | combine(cleanreads_ch)
    | STAR_align
    | map { n -> [ n.getAt(0), n.getAt(1)]}
  }

  // mapped_ch = bwamem2_mem.out // probably change this to aligned_ch
  //   | map { n -> [n.simpleName.replaceFirst("_mapped",""), n] }

  unmapped_ch = MarkIlluminaAdapters.out
    | map { n -> [n.simpleName.replaceFirst("_marked",""), n] }

  genome_ch 
    | (CreateSequenceDictionary & samtools_faidx )

  if (params.seq == "dna") {
    MergeBamAlignment_ch = unmapped_ch 
    | join(mapped_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | MergeBamAlignment_DNA
  } else if (params.seq == 'rna' ) {
    MergeBamAlignment_ch = unmapped_ch 
    | join(mapped_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | MergeBamAlignment_RNA
  }

  if(params.invariant) {
    allbambai_ch = MergeBamAlignment.out // do these need to be merged by read?
  } else {
    if(params.seq == 'dna') {
      allbai_ch = MergeBamAlignment_ch
      | map { n -> n.getAt(1)}
      | collect 
      | map { n -> [n]}

      allbambai_ch = MergeBamAlignment_ch
      | map { n -> n.getAt(0)}
      | collect
      | map { n -> [n]}
      | combine(allbai_ch)
    } else if (params.seq == 'rna') {
      MergeBamAlignment_ch
      | MarkDuplicates
      | combine(genome_ch)
      | combine(samtools_faidx.out)
      | combine(CreateSequenceDictionary.out)
      | SplitNCigarReads

      allbai_ch = SplitNCigarReads.out
      | map { n -> n.getAt(2)}
      | collect 
      | map { n -> [n]}

      allbambai_ch = SplitNCigarReads.out
      | map { n -> n.getAt(1)}
      | collect
      | map { n -> [n]}
      | combine(allbai_ch)
    }
  }

  // == Run Gatk Haplotype by interval window
  part1_ch = samtools_faidx.out
    | bedtools_makewindows
    | splitText(){it.trim()}
    | combine(allbambai_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)
  if(params.seq == 'rna'){
    part1_reads = MarkDuplicates.out
    | combine(samtools_faidx.out)
  } else {
    part1_reads = part1_ch
  }
  // If invariant, stop at output.vcf (all sites). If not, only keep SNPs with SNP filtering.
  if(params.invariant){

    part2_ch = part1_ch
      | gatk_HaplotypeCaller_invariant
      | collect
      | map { n -> [n] }
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | CombineGVCFs
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | GenotypeGVCFs

  }else{

    part2_ch = part1_ch
      | gatk_HaplotypeCaller
      | collect
      | merge_vcf
      | vcftools_snp_only
      | combine(CreateSequenceDictionary.out)
      | SortVcf
      | calc_DPvalue

    // == Filter resulting SNPs
    SortVcf.out
      | combine(calc_DPvalue.out.map{n-> n.replaceAll("\n","")})
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | VariantFiltration
      | keep_only_pass
  }

}
