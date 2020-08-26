#! /usr/bin/env nextflow

/* biocontainers are not nextflow compatable!?
fastqc_container = 'biocontainers/fastqc'
*/

fastqc_container = 'pegi3s/fastqc'

process fastqc_help {
  label 'fastqc'

  container = "$fastqc_container"
  
  output: path 'fastqc_out.txt'

  """
  fastqc --help > fastqc_out.txt
  """
}

process fastqc {
    tag "$name"
    label 'process_medium'
    label 'fastqc'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

//    when:
//    !params.skipQC && !params.skipFastQC

    input:
//    tuple val(name), file(reads) //from raw_reads_fastqc
    path(reads)

    output:
    file "*_fastqc.{zip,html}" //into fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}
