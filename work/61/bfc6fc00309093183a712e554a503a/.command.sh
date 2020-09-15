#!/bin/bash -ue
echo "I can see b73_chr1_150000001-151000000.fasta" > gatk0_index.out
bash gatk0_index.sh "b73_chr1_150000001-151000000.fasta" "indexed_genome"
