#!/bin/bash
#SBATCH --job-name=drms_perm
#SBATCH --output=drms_perm_%j.out
#SBATCH --error=drms_perm_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G

# Slurm launcher for DRMS_StTAWaveform_Permutation.m
#
# Before submitting:
#   1. Place this file and DRMS_StTAWaveform_Permutation.m in, or submit
#      from, the directory containing the MATLAB analysis and input MAT files.
#   2. Make sure the two input filenames in the MATLAB script match your files.
#   3. Change the MATLAB module name below if required by your cluster.
#   4. The MATLAB script currently uses numWorkers = 40, so this job requests
#      40 CPUs. If you change numWorkers, change --cpus-per-task to match.
#   5. Add cluster-specific #SBATCH --partition and/or --account lines above
#      if your system requires them.
#
# Submit with:
#   sbatch run_DRMS_StTAWaveform_Permutation.sbatch
#
# Check job status with:
#   squeue -u "$USER"
#
# Cancel a job with:
#   scancel JOB_ID

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

# Replace with the MATLAB module/version available on your cluster.
module load matlab

matlab -batch "DRMS_StTAWaveform_Permutation"
