#!/bin/bash
#SBATCH --job-name=Mito_Pipeline_Manager # Job name in the queue
#SBATCH --cpus-per-task=2               # Request 2 CPU cores
#SBATCH --mem=8G                          # Request 8 GB of memory
#SBATCH --time=24:00:00                   # Max runtime of 24 hours
#SBATCH --output=mito_pipeline_manager_%j.log # Log file name (%j = job ID)

# Exit immediately if any command fails
set -e

# ==============================================================================
# --- User Configuration Area (Please modify for your environment) ---
# ==============================================================================

# 1. File and Directory Paths
MAIN_NF_SCRIPT="/home/lt692/project_pi_njl27/lt692/human_mito_calling-main/main.nf"       # Path to the main Nextflow script
NEXTFLOW_CONFIG="/home/lt692/project_pi_njl27/lt692/human_mito_calling-main/nextflow.config" # Path to the Nextflow config file
MASTER_SAMPLE_LIST="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/test_sample.tsv" # TSV file containing all samples and URLs
HAIL_SCRIPT_DIR="/home/lt692/project_pi_njl27/lt692/human_mito_calling-main/hail"                 # Directory containing all hail/.py scripts
OUTPUT_DIR="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/results"    # Final results output directory (used as a fallback)

# --- CRITICAL: Define a file to store the shared, unique sample list ---
# This file is written by the Master (Stage 1) and read by the Workers (Stage 1)
UNIQUE_SAMPLE_FILE="unique_samples_for_job_array.list"

# 2. Slurm Job Array Configuration
CONCURRENT_SAMPLES=3 # Max number of samples (array tasks) to run simultaneously

# ==============================================================================
# -----------------        Script main logic       --------------------------
# This script operates in 3 modes:
# 1. Master Mode (Default): Run on login node. Creates the unique sample list and
#                          submits Stage 1 (array) and Stage 2 (finalize) jobs.
# 2. Worker Mode (if $SLURM_ARRAY_TASK_ID is set): Runs on compute nodes as an
#                          array task. Processes one unique sample.
# 3. Finalizer Mode (if $1 == --finalize): Runs on a compute node after all
#                          workers are done. Merges all results using Hail.
# ==============================================================================

# --- Mode 3: Finalizer Mode (Merging) ---
# This block executes if the script is called with the '--finalize' argument
# It will only be run by Slurm after the STAGE 1 job array (all workers) succeeds.
if [ "$1" == "--finalize" ]; then
    echo "========================================================"
    echo "--- STAGE 2: FINALIZER MODE - Merging all sample outputs"
    echo "========================================================"

    # Add the Hail script directory to the Python path
    # This ensures the 'import' statements in the .py scripts work
    export PYTHONPATH="${HAIL_SCRIPT_DIR}:${PYTHONPATH}"
    echo "[*] Added ${HAIL_SCRIPT_DIR} to PYTHONPATH."

    # Attempt to automatically find the Nextflow 'outdir' from the config file
    RESULTS_DIR=$(grep -oP "outdir\s*=\s*'\K[^']+" "$NEXTFLOW_CONFIG")
    if [ -z "$RESULTS_DIR" ]; then
        # Fallback to the user-defined OUTPUT_DIR if not found in the config
        echo "[!] WARN: Could not determine 'outdir' from nextflow.config. Using user-defined OUTPUT_DIR."
        RESULTS_DIR="$OUTPUT_DIR"
    fi
    echo "[*] Using results directory: $RESULTS_DIR"

    # --- Hail Step 1: Create coverage matrix ---
    echo "[*] Step 2.1: Creating Hail matrix table for coverage..."
    # Find all coverage files from Stage 1 and save their paths to a list
    find "$RESULTS_DIR" -name "*.per_base_coverage.tsv" > coverage_files.list
    # Run the Hail script to merge all files into a single MatrixTable (MT)
    python "${HAIL_SCRIPT_DIR}/annotate_coverage.py" --input-list coverage_files.list --output-mt coverage.mt
    echo "[+] Coverage matrix table created at coverage.mt"

    # --- Hail Step 2: Merge VCFs ---
    echo "[*] Step 2.2: Merging and filtering single sample VCF files..."
    # Find all final VCF files from Stage 1
    find "$RESULTS_DIR" -name "*.final.vcf.gz" > vcf_files.list
    # Run Hail script to combine VCFs, using the coverage.mt for filtering
    python "${HAIL_SCRIPT_DIR}/combine_vcfs.py" --coverage-mt coverage.mt --input-list vcf_files.list --output-file-name combined_variants
    echo "[+] Combined VCF and MatrixTable created."

    # --- Hail Step 3: Final Annotation ---
    echo "[*] Step 2.3: Adding additional annotations..."
    python "${HAIL_SCRIPT_DIR}/add_annotations.py" --input-mt combined_variants.mt --output-mt final_annotated_results.mt
    echo "[+] Final annotation complete."

    echo "[SUCCESS] Stage 2 (Finalizer) finished successfully."
    exit 0

# --- Mode 2: Worker Mode (Compute Node) ---
# This block executes if the script is running as part of a Slurm Job Array
# (i.e., $SLURM_ARRAY_TASK_ID is set by Slurm)
elif [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    echo "================================================================="
    echo "--- STAGE 1: WORKER MODE on Compute Node (Task ${SLURM_ARRAY_TASK_ID}) ---"
    echo "================================================================="
    
    # Critical check: Ensure the sample list file created by the Master exists
    if [ ! -f "$UNIQUE_SAMPLE_FILE" ]; then
        echo "[!] ERROR: Cannot find the unique sample list file: $UNIQUE_SAMPLE_FILE"
        echo "[!] This file should have been created by the Master Mode."
        exit 1
    fi

    # --- Key logic: Get the unique SAMPLE_ID for this specific array task ---
    # Slurm array tasks are 0-indexed, but file lines are 1-indexed
    TASK_INDEX=$((SLURM_ARRAY_TASK_ID + 1))
    # Get the sample ID by reading the Nth line (where N=TASK_INDEX) from the unique list
    SAMPLE_ID=$(awk -v line=$TASK_INDEX 'NR==line {print; exit}' "$UNIQUE_SAMPLE_FILE" | tr -d '\r') # tr removes potential carriage returns

    if [ -z "$SAMPLE_ID" ]; then
        echo "[!] No valid SAMPLE_ID found in $UNIQUE_SAMPLE_FILE for task ${SLURM_ARRAY_TASK_ID} (Line $TASK_INDEX)"
        exit 1
    fi
    echo "[*] Processing unique sample: ${SAMPLE_ID}"

    # 1. Define and enter a unique 'launch directory' for this task
    #    This prevents conflicts with .nextflow.log and .nextflow/ directories
    RUN_DIR="batch_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "${RUN_DIR}"
    cd "${RUN_DIR}"
    echo "[*] Changed to unique launch directory: $(pwd)"

    # 2. Define a unique 'work directory' for Nextflow's -w flag
    #    This is where Nextflow will store all intermediate (work) files
    WORK_DIR="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/nextflow_work/${RUN_DIR}"
    mkdir -p "$WORK_DIR"
    echo "[*] Using stable run/work directory: ${WORK_DIR}"

    # 3. Create a unique 'sample input file' for this task
    #    This file will be passed to Nextflow via --input
    TEMP_TSV="${WORK_DIR}/sample_${SAMPLE_ID}.tsv"
    
    # This grep is important: it finds *all* lines in the MASTER list that
    # match this *unique* sample ID. This correctly handles cases where
    # one sample has multiple input files (e.g., multiple lanes).
    grep "^${SAMPLE_ID}\s" "$MASTER_SAMPLE_LIST" > "$TEMP_TSV"
    
    echo "[*] Created/Updated stable input file: ${TEMP_TSV}"
    echo "[*] Input file contents:"
    cat "${TEMP_TSV}"

    # 4. Load required software modules
    module load Nextflow/24.10.2
    module load miniconda/24.11.3
    module load VEP/113.3-GCC-13.3.0

    # 5. Run the Nextflow pipeline
    nextflow run "$MAIN_NF_SCRIPT" \
        -c "$NEXTFLOW_CONFIG" \
        -profile cluster \
        -resume \
        -w ${WORK_DIR} \
        --input "$TEMP_TSV" 

    # Capture Nextflow's exit code
    NF_EXIT=$?

    if [ $NF_EXIT -eq 0 ]; then
        echo "Batch ${SLURM_ARRAY_TASK_ID} completed successfully for sample ${SAMPLE_ID}."
        echo "Cleaning up Nextflow work directory: ${WORK_DIR}"
        # Cleanup is commented out. Uncomment if you want to save space.
        #rm -rf "${WORK_DIR}"
    else
        echo "Batch ${SLURM_ARRAY_TASK_ID} failed (exit code $NF_EXIT) for sample ${SAMPLE_ID}."
        # If it failed, retain the work directory for debugging
        echo "Retaining work directory for next run/debugging: ${WORK_DIR}"
    fi

    echo "--- Finished Job Array Task ${SLURM_ARRAY_TASK_ID} ---"

# --- Mode 1: Master Mode (Login Node) ---
# This is the default block that executes when you first run the script
# Its only job is to prepare and submit the other jobs to Slurm.
else
    echo "========================================================"
    echo "--- STAGE 1: MASTER MODE on Login Node - Submitting jobs"
    echo "========================================================"
    
    # --- Key logic: Create the job array based on *unique* sample IDs ---

    echo "[*] Creating unique sample list from $MASTER_SAMPLE_LIST..."
    # 1. Get the first column (sample ID) from the master TSV.
    # 2. Remove the header row (case-insensitive).
    # 3. Sort and keep only unique sample IDs.
    # 4. Save this list to the shared file for the workers.
    cut -f1 "$MASTER_SAMPLE_LIST" | grep -v -i "^sample" | sort -u > "$UNIQUE_SAMPLE_FILE"

    # Count how many unique samples there are
    NUM_SAMPLES=$(wc -l < "$UNIQUE_SAMPLE_FILE")
    
    if [ "$NUM_SAMPLES" -eq 0 ]; then
        echo "[!] ERROR: Master sample list contains no valid sample IDs (after excluding header)."
        exit 1
    fi
    
    # Slurm array indices are 0-based
    ARRAY_INDEX=$((NUM_SAMPLES - 1)) 
    
    echo "[*] Found ${NUM_SAMPLES} unique samples (saved to $UNIQUE_SAMPLE_FILE)."
    echo "[*] Submitting job array (0-${ARRAY_INDEX}) to run ${CONCURRENT_SAMPLES} samples concurrently..."

    # Submit Stage 1 (the Worker jobs) as a job array
    # --parsable: Makes sbatch output *only* the job ID
    # --array=...%...: Runs tasks 0 to ARRAY_INDEX, with CONCURRENT_SAMPLES running at a time
    # "$0": The script to run is *this script itself* (which will trigger Mode 2)
    STAGE1_JOB_ID=$(sbatch --parsable --array=0-${ARRAY_INDEX}%${CONCURRENT_SAMPLES} "$0")
    echo "[*] Stage 1 (Nextflow runs) submitted to Slurm with Job Array ID: ${STAGE1_JOB_ID}"

    # Submit Stage 2 (the Finalizer job)
    # --dependency=afterok:...: CRITICAL! This job will only start *after*
    #                           the *entire* STAGE1_JOB_ID array completes successfully.
    # "$0" --finalize: Call *this script itself* again, but pass the '--finalize'
    #                  flag to trigger Mode 3.
    STAGE2_JOB_ID=$(sbatch --parsable --dependency=afterok:${STAGE1_JOB_ID} "$0" --finalize)
    echo "[*] Stage 2 (Hail merge) submitted to Slurm with Job ID: ${STAGE2_JOB_ID}"
    echo "[*] This job will start only after the entire job array ${STAGE1_JOB_ID} completes successfully."

    echo "[SUCCESS] All jobs submitted. Monitor with 'squeue -u \$USER'."
fi