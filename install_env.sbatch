#!/bin/bash

#SBATCH -J conda_install        # Job name, for easy identification in the Slurm queue (squeue)
#SBATCH -o conda_install_%j.out # Standard output log file (where 'echo' output goes)
#SBATCH -e conda_install_%j.err # Standard error log file (where errors go)
#SBATCH -c 4                    # Request 4 CPU cores. Mamba can use multiple threads for solving/downloading.
#SBATCH --mem=16G               # [Key] Request 16GB of memory. Complex conda resolves can be memory-intensive.
#SBATCH -t 1:00:00              # Set a 1-hour time limit. The job will be terminated if it runs longer.

echo "========== Job Started =========="
date # Print the start time to the output log

# --- Conda/Mamba Initialization ---
# This is critical for non-interactive shells (like Slurm jobs) to find the 'conda' command.
# This line may not be needed if 'conda init' was run, but it is safer to include it.
# Adjust the path if your miniconda3/anaconda3 is installed elsewhere.
source ~/miniconda3/etc/profile.d/conda.sh

# --- Environment Creation ---
# Execute the mamba create command to build the environment.
# 'mamba' is used as a much faster drop-in replacement for 'conda'.
#
# -n bio_pipeline_env: Creates a new environment named 'bio_pipeline_env'.
# -c bioconda:         Adds the 'bioconda' channel, which is essential for bioinformatics tools.
# -c conda-forge:      Adds the 'conda-forge' channel, a source for many dependencies.
# --yes:               Automatically answers 'yes' to the installation prompt.
#
# Specific package versions are "pinned" for reproducibility.
echo "Starting mamba create... This may take several minutes."
mamba create -n bio_pipeline_env -c bioconda -c conda-forge --yes \
  "bwa=0.7.18" \
  "samtools=1.21" \
  "picard-slim=3.1.1" \
  "gatk4=4.5.0.0" \
  "python=3.11" \
  "r-base=4.3" \
  "openjdk=17" # GATK and Picard require a specific Java (openjdk) version.

echo "========== Job Finished =========="
date # Print the end time to the output log