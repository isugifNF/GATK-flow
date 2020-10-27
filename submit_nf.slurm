#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=NextFlow_GATK
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
# --mail-user=username@email.com
# --mail-type=begin
# --mail-type=end
# --account=isu_gif_vrsc

set -e
set -u

start=`date +%s`

# === Load Modules here
# module load nextflow                         # If on Ceres HPC
module load gcc/7.3.0-xegsmw4 nextflow         # If on Nova/Condo HPC
module load singularity
NEXTFLOW=NextFlow
# NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow  # If on Atlas HPC
# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Get input size and module versions
# echo "started NextFlow.slurm: " `date` " seconds" >> LOGGER.txt
# module list >> LOGGER.txt
# ls -ltr ${INPUT}    # <= list any input files here
# === Main Program
#${NEXTFLOW} run main.nf -profile slurm,singularity -resume

# Method 1: Use --reads
${NEXTFLOW} run main.nf \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm,singularity \
  -resume

# Method 2: Use --reads_file
# ${NEXTFLOW} run main.nf \
#   --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
#  --reads_file read-path.txt \
#  -profile slurm,singularity \
#  -resume

end=`date +%s`

# === Log msgs and resource use
scontrol show job ${SLURM_JOB_ID}
echo "ran submit_nf.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
