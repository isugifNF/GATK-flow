#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --job-name=NF_GATK
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
# --mail-user=username@email.com
# --mail-type=begin
# --mail-type=end
# --account=project_name    # <= put HPC account name here, required on Atlas

set -e
set -u

start=`date +%s`

# === Load Modules here
# ==  Atlas HPC (will need a local install of nextflow)
# module load singularity
# NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow

# == Ceres HPC
# module load nextflow
# NEXTFLOW=nextflow

# == Nova HPC
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Main Program

# Method 1: Use --reads
${NEXTFLOW} run main.nf \
  --genome "test-data/ref/b73_chr1_150000001-151000000.fasta" \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  --queueSize 25 \
  -profile slurm,singularity \
  -resume
  #--account isu_gif_vrsc       #<= add this to Atlas

# Method 2: Use --reads_file
# ${NEXTFLOW} run main.nf \
#   --genome "test-data/ref/b73_chr1_150000001-151000000.fasta" \
#   --reads_file read-path.txt \
#   --queueSize 50 \
#   --account isu_gif_vrsc \
#   --outdir "GATK_Results" \
#   -profile slurm,singularity \
#   -resume

end=`date +%s`

# === Log msgs and resource use
scontrol show job ${SLURM_JOB_ID}
echo "ran submit_nf.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
