#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=NextFlow
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=jenchang@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e
set -u

start=`date +%s`

# === Load Modules here
# module load nextflow                         # If on Ceres HPC
module load gcc/7.3.0-xegsmw4 nextflow       # If on Condo HPC, don't know why this seems more complicated

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Get input size and module versions
# echo "started NextFlow.slurm: " `date` " seconds" >> LOGGER.txt
# module list >> LOGGER.txt
# ls -ltr ${INPUT}    # <= list any input files here
# === Main Program
nextflow run main.nf -profile condo -resume

end=`date +%s`

# === Log msgs and resource use                          
scontrol show job ${SLURM_JOB_ID}
echo "ran submit_nf.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
