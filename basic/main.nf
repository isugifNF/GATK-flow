#!/usr/bin/env nextflow

nextflow.enable.dsl=2
println "Hello world!"

params.n = "./*.nf"

process myhi {
  input:
  path name

  output:
  stdout() 

  script:
  """
  echo "hello! $name"
  """
}


workflow {

channel.fromPath("$params.n") | myhi

}
