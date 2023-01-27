#! /usr/bin/env nextflow

nextflow.enable.dsl=2

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
  'genome','reads','reads_file','invariant',
  'singularity_img','docker_img',
  'gatk_app','bwamem2_app','samtools_app','bedtools_app','datamash_app','vcftools_app',
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


process FastqToSam {
  tag "$readname"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // [readgroup, [left.fq.gz, right.fq.gz], increment_readgroup]
  tuple val(readname), path(readpairs), val(i_readname)

  output: // increment_readgroup.bam since one readgroup can have multiple lanes
  path("*.bam")
  
  script:
  """
  #! /usr/bin/env bash
  ${gatk_app} --java-options "${java_options}" FastqToSam \
    --FASTQ ${readpairs.get(0)} \
    --FASTQ2 ${readpairs.get(1)} \
    --OUTPUT ${i_readname}.bam \
    --READ_GROUP_NAME ${readname} \
    --SAMPLE_NAME ${readname}_name \
    --LIBRARY_NAME ${readname}_lib \
    --PLATFORM ILLUMINA \
    --SEQUENCING_CENTER ISU \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${i_readname}.bam
  """  
}

process MarkIlluminaAdapters {
  tag "${bam.fileName}"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // reads.bam
  path(bam)

  output: // reads_marked.bam
  path "${bam.simpleName}_marked.bam"

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" MarkIlluminaAdapters \
    --INPUT $bam \
    --OUTPUT ${bam.simpleName}_marked.bam \
    --METRICS ${bam.simpleName}_marked_metrics.txt \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_marked.bam
  touch ${bam.simpleName}_marked_metrics.txt
  """
}

process SamToFastq {
  tag "${bam.fileName}"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // reads.bam
  path(bam)

  output: // reads_interleaved.fq
  tuple val("${bam.simpleName}"), path("${bam.simpleName}_newR1.fq"), path("${bam.simpleName}_newR2.fq")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" SamToFastq \
    --INPUT $bam \
    --FASTQ ${bam.simpleName}_newR1.fq \
    --SECOND_END_FASTQ ${bam.simpleName}_newR2.fq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INCLUDE_NON_PF_READS true \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_newR1.fq
  touch ${bam.simpleName}_newR2.fq
  """
}
// INTERLEAVE=true
// USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

process bwamem2_index {
  tag "${genome_fasta.simpleName}"
  label 'bwamem'
  publishDir "${params.outdir}/02_MapReads"

  input:
  path(genome_fasta)

  output: // [genome.fasta, [genome_index files]]
  tuple path("$genome_fasta"), path("${genome_fasta}*")

  script:
  """
  #! /usr/bin/env bash
  $bwamem2_app index $genome_fasta
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta}.variousindexfiles
  """
}

process bwamem2_mem {
  tag "$readname"
  label 'bwamem'
  publishDir "${params.outdir}/02_MapReads"

  input:
  tuple path(genome_fasta), path(genome_index), val(readname), path(readpairs)

  output: // reads_mapped_2_genome.bam
  path("${readname}_mapped.bam")

  script:
  """
  #! /usr/bin/env bash
  PROC1=\$((`nproc` * 3/4))
  $bwamem2_app mem -t \${PROC1} ${genome_fasta} ${readpairs} |\
     $samtools_app view --threads 1 -bS - > ${readname}_mapped.bam
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${readname}_mapped.bam
  """
}

process CreateSequenceDictionary {
  tag "${genome_fasta.simpleName}"
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"

  input:
  path(genome_fasta)

  output:
  path("${genome_fasta.simpleName}.dict")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" CreateSequenceDictionary \
    -R ${genome_fasta} \
    -O ${genome_fasta.simpleName}.dict
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta.simpleName}.dict
  """
}

process MergeBamAlignment {
  tag "$i_readname"
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"

  input:  // [readgroup, unmapped reads, mapped reads]
  tuple val(i_readname), path(read_unmapped), path(read_mapped), path(genome_fasta), path(genome_dict)

  output: // merged bam and bai files
  tuple path("${i_readname}_merged.bam"), path("${i_readname}_merged.bai")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" MergeBamAlignment \
  --REFERENCE_SEQUENCE $genome_fasta \
  --UNMAPPED_BAM ${read_unmapped} \
  --ALIGNED_BAM ${read_mapped} \
  --OUTPUT ${i_readname}_merged.bam \
  --CREATE_INDEX true \
  --ADD_MATE_CIGAR true \
  --CLIP_ADAPTERS false \
  --CLIP_OVERLAPPING_READS true \
  --INCLUDE_SECONDARY_ALIGNMENTS true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --ATTRIBUTES_TO_RETAIN XS \
  --USE_JDK_DEFLATER true \
  --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${i_readname}_merged.bam
  touch ${i_readname}_merged.bai
  """
}

process samtools_faidx {
  tag "${genome_fasta.simpleName}"
  label 'samtools'
  publishDir "${params.outdir}/03_PrepGATK"

  input:
  path(genome_fasta)

  output:
  path("${genome_fasta}.fai")

  script:
  """
  #! /usr/bin/env bash
  $samtools_app faidx $genome_fasta
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta}.fai
  """
}

process bedtools_makewindows {
  tag "${genome_fai.simpleName}"
  label 'bedtools'
  publishDir "${params.outdir}/03_PrepGATK"

  input:  // genome.fai
  path(genome_fai)

  output: // genome_coords.bed
  path("*_coords.bed")

  script:
  """
  #! /usr/bin/env bash
  cat ${genome_fai} | awk -F'\t' '{print \$1"\t"\$2}' > genome_length.txt
  $bedtools_app makewindows -w "$params.window" -g genome_length.txt |\
    awk '{print \$1"\t"\$2+1"\t"\$3}' |\
    sed \$'s/\t/:/1' |\
    sed \$'s/\t/-/g' > ${genome_fai.simpleName}_coords.bed
  """

  stub:
  """
  #! /usr/bin/env bash
  echo "0-1000" > ${genome_fai.simpleName}_coords.bed
  echo "1000-2000" >> ${genome_fai.simpleName}_coords.bed
  echo "2000-3000" >> ${genome_fai.simpleName}_coords.bed
  """
}
//  ${fai.simpleName}_coords.bed


process gatk_HaplotypeCaller {
  tag "$window"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple val(window), path(bam), path(bai), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  BAMFILES=`echo $bam | sed 's/ / -I /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options "${java_options}" HaplotypeCaller \
    -R $genome_fasta \
    -I \$BAMFILES \
    -L $window \
    --output ${window.replace(':','_')}.vcf
  """

  stub:
  """
  touch ${window.replace(':','_')}.vcf
  """
}

process gatk_HaplotypeCaller_invariant {
  tag "$window:${bam.simpleName}"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple val(window), path(bam), path(bai), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  BAMFILES=`echo $bam | sed 's/ / -I /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options "${java_options}" HaplotypeCaller \
    -ERC BP_RESOLUTION \
    -R $genome_fasta \
    -I \$BAMFILES \
    -L $window \
    --output ${bam.simpleName}_${window.replace(':','_')}.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_${window.replace(':','_')}.vcf
  """
}
// --include-invariant -ERC GVCF

// Consolidate GVCFs
process CombineGVCFs {
  tag "ConsolidateGVCFs"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple path(gvcf), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  GVCFFILES=`echo $gvcf | sed 's/ / --variant /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options "${java_options}" CombineGVCFs \
    -R $genome_fasta \
    --variant \$GVCFFILES \
    --output all_combined.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch all_combined.vcf
  """
}

// Joint-Call Cohorts
process GenotypeGVCFs {
  tag "JointCallCohorts"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple path(all_combined_gvcf), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" GenotypeGVCFs \
    -R $genome_fasta \
    -V $all_combined_gvcf \
    --output output.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch output.vcf
  """
}

process merge_vcf {
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // multiple SNP vcf files
  path(vcfs)

  output: // merged into one vcf file
  path("first-round_merged.vcf")

  script:
  """
  #! /usr/bin/env bash
  cat ${vcfs.get(0)} | grep "^#" > first-round_merged.vcf
  cat ${vcfs} | grep -v "^#" >> first-round_merged.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch first-round_merged.vcf
  """
}

process vcftools_snp_only {
  tag "${merged_vcf.fileName}"
  label 'vcftools'
  publishDir "${params.outdir}/05_FilterSNPs", mode: 'copy'

  input:  // merged SNP vcf file
  path(merged_vcf)

  output: // vcf file only containing SNPs
  path("${merged_vcf.simpleName}_snps-only.*")

  script:
  """
  #! /usr/bin/env bash
  $vcftools_app \
    --vcf $merged_vcf \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --out ${merged_vcf.simpleName}_snps-only
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${merged_vcf.simpleName}_snps-only.outputfiles
  """
}

process SortVcf {
  tag "$vcf.fileName"
  label 'gatk'
  publishDir "$params.outdir/05_FilterSNPs", mode: 'copy'

  input:  // [SNP.vcf, genome.dict]
  tuple path(vcf), path(dict)

  output: // sorted SNP.vcf
  path ("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" SortVcf \
  --INPUT $vcf \
  --SEQUENCE_DICTIONARY $dict \
  --CREATE_INDEX true \
  --OUTPUT ${vcf.simpleName}_sorted.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${vcf.simpleName}_sorted.vcf
  """
}
// java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar

process calc_DPvalue {
  tag "$sorted_vcf.fileName"
  label 'datamash'

  input:  // sorted SNP vcf
  path(sorted_vcf)

  output: // DP value (number) to filter vcf in downstream process
  stdout()

  script:
  """
  #! /usr/bin/env bash
  grep -v "^#" $sorted_vcf | cut -f 8 | grep -oe ";DP=.*" | cut -f 2 -d ';' | cut -f 2 -d "=" > dp.txt
  cat dp.txt | $datamash_app mean 1 sstdev 1 > dp.stats
  cat dp.stats | awk '{print \$1+5*\$2}'
  """

  stub:
  """
  #! /usr/bin/env bash
  # idk what default number should be here
  echo "75"
  """
}

process VariantFiltration {
  tag "$sorted_snp_vcf.fileName"
  publishDir "$params.outdir/05_FilterSNPs", mode: 'copy'

  input:  // [sorted snp vcf, DP filter, genome files ... ]
  tuple path(sorted_snp_vcf), val(dp), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // filtered to identified SNP variants
  path("${sorted_snp_vcf.simpleName}.marked.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" VariantFiltration \
    --reference $genome_fasta \
    --sequence-dictionary $genome_dict \
    --variant $sorted_snp_vcf \
    --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > $dp\" \
    --filter-name "FAIL" \
    --output ${sorted_snp_vcf.simpleName}.marked.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${sorted_snp_vcf.simpleName}.marked.vcf
  """
}
// java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar

process keep_only_pass {
  tag "${snp_marked_vcf.fileName}"
  publishDir "$params.outdir/05_FilterSNPs", mode: 'copy'

  input:
  path(snp_marked_vcf)

  output:
  path("${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf")

  script:
  """
  #! /usr/bin/env bash
  cat $snp_marked_vcf | grep "^#" > ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  cat $snp_marked_vcf | grep -v "^#" | awk '\$7=="PASS"' >> ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  """
}

workflow {
  // == Read in genome and reads channels
  if(params.genome) {
    genome_ch = channel.fromPath(params.genome, checkIfExists:true)
      | view {file -> "Genome file : $file "}
  } else {
    exit 1, "[Missing File(s) Error] Maize_WGS_Build requires a reference '--genome [GENOME.fasta]' \n"
  }

  if (params.reads) {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
      | view {files -> "Read files : $files "}
  } else if (params.reads_file) {
    reads_ch = channel.fromPath(params.reads_file, checkIfExists:true)
      | splitCsv(sep:'\t')
      | map { n -> [ n.get(0), [n.get(1), n.get(2)]] }
      | view {files -> "Read files : $files "}
  } else {
    exit 1, "[Missing File(s) Error] Maize_WGS_Build requires either paired-end read files as a glob '--reads [*_{r1,r2}.fq.gz]' or as a tab-delimited text file '--reads_file [READS_FILE.txt]'\n"
  }

  // == Since one sample may be run on multiple lanes
  i = 1

  // == Prepare mapped and unmapped read files
  cleanreads_ch = reads_ch
    | map { n -> [n.get(0), n.get(1), "${i++}_"+n.get(0)] }
    | FastqToSam
    | MarkIlluminaAdapters
    | SamToFastq
    | map { n -> [ n.get(0).replaceFirst("_marked",""), [ n.get(1), n.get(2)] ] }

  genome_ch
    | bwamem2_index
    | combine(cleanreads_ch)
    | bwamem2_mem

  mapped_ch = bwamem2_mem.out
    | map { n -> [n.simpleName.replaceFirst("_mapped",""), n] }

  unmapped_ch = MarkIlluminaAdapters.out
    | map { n -> [n.simpleName.replaceFirst("_marked",""), n] }

  genome_ch 
    | (CreateSequenceDictionary & samtools_faidx )
  unmapped_ch 
    | join(mapped_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | MergeBamAlignment

  if(params.invariant) {
    allbambai_ch = MergeBamAlignment.out // do these need to be merged by read?
  } else {
    allbai_ch = MergeBamAlignment.out
      | map { n -> n.get(1)}
      | collect | map { n -> [n]}

    allbambai_ch = MergeBamAlignment.out
      | map { n -> n.get(0)}
      | collect
      | map { n -> [n]}
      | combine(allbai_ch)
  }

  // == Run Gatk Haplotype by interval window
  part1_ch = samtools_faidx.out
    | bedtools_makewindows
    | splitText(){it.trim()}
    | combine(allbambai_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)

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
