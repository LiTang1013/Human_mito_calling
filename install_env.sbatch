#!/bin/bash

#SBATCH -J conda_install      # 任务名，方便识别
#SBATCH -p general            # 提交到名为 general 的分区 (请根据您集群的实际情况修改)
#SBATCH -o conda_install_%j.out # 标准输出日志
#SBATCH -e conda_install_%j.err # 错误输出日志
#SBATCH -c 4                  # 申请4个CPU核心
#SBATCH --mem=16G             # 【关键】申请16GB内存，这应该远超安装所需
#SBATCH -t 1:00:00            # 任务最长运行1小时

echo "========== Job Started =========="
date

# 确保Conda已经初始化
# 这行根据您.bashrc的配置可能需要，也可能不需要，加上更保险
source ~/miniconda3/etc/profile.d/conda.sh

# 执行我们之前定义的安装命令
mamba create -n bio_pipeline_env -c bioconda -c conda-forge --yes \
  "bwa=0.7.18" \
  "samtools=1.21" \
  "picard-slim=3.1.1" \
  "gatk4=4.5.0.0" \
  "python=3.11" \
  "r-base=4.3" \
  "openjdk=17"

echo "========== Job Finished =========="
date