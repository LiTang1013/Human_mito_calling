#!/bin/bash
#SBATCH --job-name=Mito_Pipeline_Manager
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00 
#SBATCH --output=mito_pipeline_manager_%j.log

set -e

# ==============================================================================
# --- 用户配置区 (请根据您的环境修改) ---
# =-============================================================================

# 1. 文件与目录路径
MAIN_NF_SCRIPT="/home/lt692/project_pi_njl27/lt692/human_mito_calling-main/main.nf"
NEXTFLOW_CONFIG="/home/lt692/project_pi_njl27/lt692/human_mito_calling-main/nextflow.config"
MASTER_SAMPLE_LIST="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/test_sample.tsv" # 包含所有样本和URL的TSV文件
HAIL_SCRIPT_DIR="/home/lt692/project_pi_njl27/lt692/human_mito_calling-main/hail"                 # 包含所有hail/.py脚本的目录
OUTPUT_DIR="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/results"    # 最终结果输出目录

# --- 关键：定义一个文件，用于存储 Master 和 Worker 之间共享的唯一样本列表 ---
UNIQUE_SAMPLE_FILE="unique_samples_for_job_array.list"

# 2. Slurm 作业阵列配置
CONCURRENT_SAMPLES=3

# ==============================================================================
# --- 脚本主逻辑 (通常无需修改) ---
# ==============================================================================

# --- 模式三：最终合并阶段 (Finalizer Mode) ---
if [ "$1" == "--finalize" ]; then
    echo "========================================================"
    echo "--- STAGE 2: FINALIZER MODE - Merging all sample outputs"
    echo "========================================================"

    export PYTHONPATH="${HAIL_SCRIPT_DIR}:${PYTHONPATH}"
    echo "[*] Added ${HAIL_SCRIPT_DIR} to PYTHONPATH."

    RESULTS_DIR=$(grep -oP "outdir\s*=\s*'\K[^']+" "$NEXTFLOW_CONFIG")
    if [ -z "$RESULTS_DIR" ]; then
        echo "[!] WARN: Could not determine 'outdir' from nextflow.config. Using user-defined OUTPUT_DIR."
        RESULTS_DIR="$OUTPUT_DIR"
    fi
    echo "[*] Using results directory: $RESULTS_DIR"

    # --- Hail 步骤 1: 创建覆盖度矩阵 ---
    echo "[*] Step 2.1: Creating Hail matrix table for coverage..."
    find "$RESULTS_DIR" -name "*.per_base_coverage.tsv" > coverage_files.list
    python "${HAIL_SCRIPT_DIR}/annotate_coverage.py" --input-list coverage_files.list --output-mt coverage.mt
    echo "[+] Coverage matrix table created at coverage.mt"

    # --- Hail 步骤 2: 合并VCF ---
    echo "[*] Step 2.2: Merging and filtering single sample VCF files..."
    find "$RESULTS_DIR" -name "*.final.vcf.gz" > vcf_files.list
    python "${HAIL_SCRIPT_DIR}/combine_vcfs.py" --coverage-mt coverage.mt --input-list vcf_files.list --output-file-name combined_variants
    echo "[+] Combined VCF and MatrixTable created."

    # --- Hail 步骤 3: 最终注释 ---
    echo "[*] Step 2.3: Adding additional annotations..."
    python "${HAIL_SCRIPT_DIR}/add_annotations.py" --input-mt combined_variants.mt --output-mt final_annotated_results.mt
    echo "[+] Final annotation complete."

    echo "[SUCCESS] Stage 2 finished successfully."
    exit 0

# --- 模式二：计算节点工作模式 (Worker Mode) ---
elif [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    echo "================================================================="
    echo "--- STAGE 1: WORKER MODE on Compute Node (Task ${SLURM_ARRAY_TASK_ID}) ---"
    echo "================================================================="
    
    # 检查由 Master Mode 创建的唯一样本列表是否存在
    if [ ! -f "$UNIQUE_SAMPLE_FILE" ]; then
        echo "[!] ERROR: Cannot find the unique sample list file: $UNIQUE_SAMPLE_FILE"
        echo "[!] This file should have been created by the Master Mode."
        exit 1
    fi

    # --- 关键修复：从 UNIQUE_SAMPLE_FILE 中按行号提取唯一的 SAMPLE_ID ---
    # (不再使用 awk 从 MASTER_SAMPLE_LIST 中提取)
    TASK_INDEX=$((SLURM_ARRAY_TASK_ID + 1))
    SAMPLE_ID=$(awk -v line=$TASK_INDEX 'NR==line {print; exit}' "$UNIQUE_SAMPLE_FILE" | tr -d '\r')

    if [ -z "$SAMPLE_ID" ]; then
        echo "[!] No valid SAMPLE_ID found in $UNIQUE_SAMPLE_FILE for task ${SLURM_ARRAY_TASK_ID} (Line $TASK_INDEX)"
        exit 1
    fi
    echo "[*] Processing unique sample: ${SAMPLE_ID}"

    # 1. 定义并进入唯一的“启动目录” (e.g., batch_0)
    #    这可以防止 .nextflow.log 和 .nextflow/ 目录冲突
    RUN_DIR="batch_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "${RUN_DIR}"
    cd "${RUN_DIR}"
    echo "[*] Changed to unique launch directory: $(pwd)"

    # 2. 定义唯一的“工作目录” (用于 -w, e.g., /scratch/.../batch_0)
    WORK_DIR="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/nextflow_work/${RUN_DIR}"
    mkdir -p "$WORK_DIR"
    echo "[*] Using stable run/work directory: ${WORK_DIR}"

    # 3. 创建唯一的“样本输入文件” (在 WORK_DIR 中)
    #    *** 这一步的逻辑是正确的 ***
    #    它会从 MASTER_SAMPLE_LIST 中 grep *所有* 匹配该 SAMPLE_ID 的行
    TEMP_TSV="${WORK_DIR}/sample_${SAMPLE_ID}.tsv"
    grep "^${SAMPLE_ID}\s" "$MASTER_SAMPLE_LIST" > "$TEMP_TSV"
    echo "[*] Created/Updated stable input file: ${TEMP_TSV}"
    echo "[*] Input file contents:"
    cat "${TEMP_TSV}"

    # 4. 加载模块
    module load Nextflow/24.10.2
    module load miniconda/24.11.3
    module load VEP/113.3-GCC-13.3.0

    # 5. 运行 Nextflow
    nextflow run "$MAIN_NF_SCRIPT" \
        -c "$NEXTFLOW_CONFIG" \
        -profile cluster \
        -resume \
        -w ${WORK_DIR} \
        --input "$TEMP_TSV" 

    NF_EXIT=$?

    if [ $NF_EXIT -eq 0 ]; then
        echo "Batch ${SLURM_ARRAY_TASK_ID} completed successfully for sample ${SAMPLE_ID}."
        echo "Cleaning up Nextflow work directory: ${WORK_DIR}"
        #rm -rf "${WORK_DIR}"
    else
        echo "Batch ${SLURM_ARRAY_TASK_ID} failed (exit code $NF_EXIT) for sample ${SAMPLE_ID}."
        echo "Retaining work directory for next run/debugging: ${WORK_DIR}"
    fi

    echo "--- Finished Job Array Task ${SLURM_ARRAY_TASK_ID} ---"

# --- 模式一：登录节点主控模式 (Master Mode) ---
else
    echo "========================================================"
    echo "--- STAGE 1: MASTER MODE on Login Node - Submitting jobs"
    echo "========================================================"
    
    # --- 关键修复：基于唯一的样本ID创建作业阵列 ---

    echo "[*] Creating unique sample list from $MASTER_SAMPLE_LIST..."
    # 1. 从 MASTER_SAMPLE_LIST 的第一列 (cut -f1)
    # 2. 移除标题行 (grep -v -i "^sample")
    # 3. 排序并只保留唯一的 (sort -u)
    # 4. 保存到共享文件
    cut -f1 "$MASTER_SAMPLE_LIST" | grep -v -i "^sample" | sort -u > "$UNIQUE_SAMPLE_FILE"

    # 统计唯一样本的数量
    NUM_SAMPLES=$(wc -l < "$UNIQUE_SAMPLE_FILE")
    
    if [ "$NUM_SAMPLES" -eq 0 ]; then
        echo "[!] ERROR: Master sample list contains no valid sample IDs (after excluding header)."
        exit 1
    fi
    # Slurm 数组索引是从 0 开始的
    ARRAY_INDEX=$((NUM_SAMPLES - 1)) 
    
    echo "[*] Found ${NUM_SAMPLES} unique samples (saved to $UNIQUE_SAMPLE_FILE)."
    echo "[*] Submitting job array (0-${ARRAY_INDEX}) to run ${CONCURRENT_SAMPLES} samples concurrently..."

    # 提交阶段一的作业阵列，并捕获其Job ID
    STAGE1_JOB_ID=$(sbatch --parsable --array=0-${ARRAY_INDEX}%${CONCURRENT_SAMPLES} "$0")
    echo "[*] Stage 1 (Nextflow runs) submitted to Slurm with Job Array ID: ${STAGE1_JOB_ID}"

    # 提交阶段二的作业，并设置其依赖于阶段一的作业阵列
    STAGE2_JOB_ID=$(sbatch --parsable --dependency=afterok:${STAGE1_JOB_ID} "$0" --finalize)
    echo "[*] Stage 2 (Hail merge) submitted to Slurm with Job ID: ${STAGE2_JOB_ID}"
    echo "[*] This job will start only after the entire job array ${STAGE1_JOB_ID} completes successfully."

    echo "[SUCCESS] All jobs submitted. Monitor with 'squeue -u \$USER'."
fi