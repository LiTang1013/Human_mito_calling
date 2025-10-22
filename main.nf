#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ====================================================================================
//                                  PIPELINE HEADER
// ====================================================================================
log.info """
PIPELINE START
=================================
Input file:        ${params.input}
Output directory:  ${params.outdir}
Reference Genome:  ${params.ref_genome_fasta}
WDL Script:        ${params.wdl_script ?: 'N/A'}
=================================
"""

// ====================================================================================
//                                  INPUT CHANNELS
// ====================================================================================
ch_ref_genome_tuple = tuple( file(params.ref_genome_fasta), file("${params.ref_genome_fasta}.{amb,ann,bwt,pac,sa,fai}") )
ch_fasta_only = Channel.value( ch_ref_genome_tuple[0] )
ch_cromwell_conf = file("${baseDir}/cromwell.conf")
ch_hail_directory = file("${baseDir}/hail")

// --- Channel for input samples ---
ch_samples = Channel.fromPath(params.input)
    .splitCsv(header: false, sep: '\t') 
    .filter { row -> row.size() == 2 && row[0] && !row[0].trim().isEmpty() && row[1] && !row[1].trim().isEmpty() } 
    .map { row -> tuple(row[0], row[1]) } 
    .groupTuple() 
    .map { sample_id, urls_list ->
        def meta = [id: sample_id]
        def sorted_urls = urls_list.sort()
        if (sorted_urls.size() % 2 != 0) {
            error "Sample ${meta.id} has an odd number of FASTQ files after grouping. Please check input file."
        }
        def fastq_pairs = []
        for (int i = 0; i < sorted_urls.size(); i += 2) {
            fastq_pairs.add([sorted_urls[i], sorted_urls[i+1]])
        }
        return [meta, fastq_pairs]
    }
    .flatMap { meta, fastq_pairs ->
        // 将每个 FASTQ 对拆分成一个独立的任务，并添加 pair_id
        fastq_pairs.indexed().collect { i, pair -> [meta + [pair_id: i], pair[0], pair[1]] }
    }

// ====================================================================================
//                                    WORKFLOW (通道重路由)
// ====================================================================================
workflow {
    DOWNLOAD_FASTQ(ch_samples) // Output: [meta (with pair_id), fq1, fq2]

    // 1. 对每对 FASTQ 进行单独比对 (输出未排序 BAM)
    ALIGN_AND_UNSORT(DOWNLOAD_FASTQ.out, ch_ref_genome_tuple) // Output: [meta, unsorted_bam]

    // 2. 排序、转换为 CRAM、索引并清理 FASTQ
    SORT_AND_CONVERT_TO_CRAM(ALIGN_AND_UNSORT.out.unsorted_bam, ch_fasta_only) // Output: [meta, single_cram, single_crai]

    // 3. 收集该样本的所有 CRAM 文件 (用于合并)
    ch_crams_to_merge = SORT_AND_CONVERT_TO_CRAM.out.single_cram
        .map { meta, cram, crai -> [meta.id, [cram, crai]] } // 映射为 [sample_id, [cram, crai]]
        .groupTuple() // 按照 sample_id 分组: [sample_id, [[cram1, crai1], [cram2, crai2], ...]]
        .map { sample_id, cram_pair_list ->
            def meta = [id: sample_id]
            // 分别提取 cram 文件列表和 crai 文件列表
            def crams = cram_pair_list.collect { it[0] }
            def crais = cram_pair_list.collect { it[1] }
            return [meta, crams, crais] // 重新包装为 [meta, [cram1, ...], [crai1, ...]]
        }

    // 4. 合并所有 CRAM 文件并删除中间文件
    MERGE_CRAMS(ch_crams_to_merge) // Output: [meta, merged_cram, merged_crai]

    // 5. 后续流程连接到 MERGE_CRAMS 的输出
    GENERATE_CRAM_TSV(MERGE_CRAMS.out.merged_cram)
    GENERATE_WDL_JSON(GENERATE_CRAM_TSV.out.tsv)
    RUN_WDL_VARIANT_CALLING(GENERATE_WDL_JSON.out.json, ch_cromwell_conf)

    ANNOTATE_INDIVIDUAL_VCF(RUN_WDL_VARIANT_CALLING.out.wdl_results, ch_hail_directory)
}

// --- Workflow-level event handlers for reporting ---
workflow.onComplete {
    log.info "Pipeline completed successfully. Output files are in: ${params.outdir}"
}
workflow.onError {
    log.error """
    ----------------------------------------------------
    ERROR: Pipeline execution failed!
    The last error message was: ${task.errorReport}
    Process name:           ${task.name}
    Process tag:            ${task.tag}
    ----------------------------------------------------
    """
}

// ====================================================================================
//                                    PROCESSES (修改后)
// ====================================================================================

process DOWNLOAD_FASTQ {
    tag "Download for ${meta.id} (Pair ${meta.pair_id})"

    // 仍然发布到结果目录，但 ALIGN_AND_UNSORT 成功后会删除工作目录中的副本
    publishDir "${params.outdir}/${meta.id}/fastq", mode: 'copy' 

    input:
    tuple val(meta), val(url1), val(url2)
    output:
    tuple val(meta), path("${meta.id}_${meta.pair_id}_1.fastq.gz"), path("${meta.id}_${meta.pair_id}_2.fastq.gz")
    script:
    """
    #!/bin/bash
    set -e
    # Wget 下载
    wget --no-check-certificate -O "${meta.id}_${meta.pair_id}_1.fastq.gz" "${url1}"
    wget --no-check-certificate -O "${meta.id}_${meta.pair_id}_2.fastq.gz" "${url2}"
    """
}


process ALIGN_AND_UNSORT {
    tag "BWA-MEM on ${meta.id} (Pair ${meta.pair_id})"

    // 资源请求应调整以匹配 BWA-MEM 自身的需求，通常比排序低
    cpus 8 
    memory '32 GB'
    time '1 h'

    input:
    // 接收单个 FASTQ 对
    tuple val(meta), path(read1), path(read2)
    tuple path(ref_fasta), path(bwa_indices)
    
    output:
    // 输出单个未排序的 BAM 文件
    tuple val(meta), path("${meta.id}_${meta.pair_id}.unsorted.bam"), emit: unsorted_bam 
    
    script:
    // 关键：Read Group ID 必须唯一，使用 meta.id 和 pair_id 组合
    def read_group_id = "${meta.id}_${meta.pair_id}" 
    def read_group = "\'@RG\\tID:${read_group_id}\\tSM:${meta.id}\\tPL:ILLUMINA\'"
    """
    #!/bin/bash
    set -e

    # 只执行比对，输出未排序的 BAM
    bwa mem -t ${task.cpus} -R ${read_group} ${ref_fasta} ${read1} ${read2} | \\
    samtools view -@ ${task.cpus} -b -o ${meta.id}_${meta.pair_id}.unsorted.bam -
    
    # 成功生成 BAM 文件后，删除原始的 FASTQ 文件 (仅删除工作目录中的副本)
    rm ${read1} ${read2}
    """
}

process SORT_AND_CONVERT_TO_CRAM {
    tag "Sort and CRAM for ${meta.id} (Pair ${meta.pair_id})"

    // 增加资源请求以解决 OOM 问题
    cpus 16
    memory '64 GB'
    time '1 h'
    
    input:
    // 接收单个未排序的 BAM 文件
    tuple val(meta), path(unsorted_bam)
    path ref_fasta
    
    output:
    // 输出单个 CRAM 和索引
    tuple val(meta), path("${meta.id}_${meta.pair_id}.cram"), path("${meta.id}_${meta.pair_id}.cram.crai"), emit: single_cram
    
    script:
    // Groovy 变量（仅用于文件名构造，不直接用于 Bash 逻辑）
    def unsorted_bam_filename = unsorted_bam.getName()

    """
    #!/bin/bash
    set -e
    
    # 1. 定义 Bash 变量
    SORTED_BAM="${meta.id}_${meta.pair_id}.sorted.bam"
    OUTPUT_CRAM="${meta.id}_${meta.pair_id}.cram"

    # 2. 排序：使用 samtools sort
    # 注意：这里我们使用 Nextflow 传递进来的文件名 ${unsorted_bam_filename}
    samtools sort -@ ${task.cpus} -o \${SORTED_BAM} ${unsorted_bam_filename}
    
    # 3. 转换为 CRAM 格式
    samtools view -@ ${task.cpus} -T ${ref_fasta} -C -o \${OUTPUT_CRAM} \${SORTED_BAM}

    # 4. 创建 CRAM 索引
    samtools index \${OUTPUT_CRAM}

    # 5. 清理：删除中间文件 (未排序和已排序的 BAM)
    rm ${unsorted_bam_filename} \${SORTED_BAM}
    """
}

process MERGE_CRAMS {
    tag "Merge CRAMs for ${meta.id}"

    // 发布最终的 CRAM 文件
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy', pattern: "*.{cram,crai}"

    input:
    // 接收该样本的所有 CRAM 文件列表
    tuple val(meta), path(crams), path(crais)

    output:
    // 输出合并后的 CRAM 和 CRAI
    tuple val(meta), path("${meta.id}.merged.cram"), path("${meta.id}.merged.cram.crai"), emit: merged_cram

    script:
    def cram_files = crams.join(' ')
    def crai_files = crais.join(' ')

    """
    #!/bin/bash
    set -e
    
    # 1. 定义 Bash 变量，确保输出文件名正确
    OUTPUT_CRAM="${meta.id}.merged.cram"

    # 2. 使用 samtools merge 合并所有 CRAM 文件
    # 注意：-f 强制覆盖，-O CRAM 指定输出格式
    samtools merge -@ ${task.cpus} -f -O CRAM --output-fmt-option version=3.0 -o \${OUTPUT_CRAM} ${cram_files}

    # 3. 创建 CRAM 索引
    samtools index \${OUTPUT_CRAM}
    
    # 4. 清理：删除所有原始的单对 CRAM 文件和它们的索引
    rm ${cram_files} ${crai_files}
    """
}

// ... (GENERATE_CRAM_TSV, GENERATE_WDL_JSON, RUN_WDL_VARIANT_CALLING, ANNOTATE_INDIVIDUAL_VCF 进程保持不变)

process GENERATE_CRAM_TSV {
    tag "Generate CRAM TSV for ${meta.id}"

    publishDir "${params.outdir}/${meta.id}/variant_calling/inputs", mode: 'copy'
    input:
    tuple val(meta), path(cram), path(crai)
    output:
    tuple val(meta), path("${meta.id}_cram_list.tsv"), emit: tsv
    script:
    """
    cram_path=\$(readlink -f ${cram})
    crai_path=\$(readlink -f ${crai})
    echo -e "\${cram_path}\\t\${crai_path}" > ${meta.id}_cram_list.tsv
    """
}

// ... (GENERATE_WDL_JSON, RUN_WDL_VARIANT_CALLING, ANNOTATE_INDIVIDUAL_VCF 保持不变)

process GENERATE_WDL_JSON {
    tag "Generate WDL JSON for ${meta.id}"
    publishDir "${params.outdir}/${meta.id}/variant_calling/inputs", mode: 'copy'

    input:
    tuple val(meta), path(cram_tsv)

    output:
    tuple val(meta), path("${meta.id}_wdl_inputs.json"), emit: json

    script:
    // 1. In Groovy, create a JSON template with a unique placeholder.
    def json_content = new groovy.json.JsonBuilder()
    def wdl_inputs = params.wdl_inputs.collectEntries { key, value ->
        ["MitochondriaMultiSamplePipeline.${key}", value]
    }
    wdl_inputs["MitochondriaMultiSamplePipeline.inputSamplesFile"] = "___TSV_PATH_PLACEHOLDER___"
    json_content(wdl_inputs)
    def json_string_template = json_content.toPrettyString()

    """
    #!/bin/bash
    set -e

    # 2. In Bash, get the absolute path of the TSV file.
    TSV_PATH=\$(readlink -f ${cram_tsv})

    # 3. Use 'printf' to create the JSON template string in a shell variable.
    JSON_TEMPLATE=\$(printf '%s' '${json_string_template}')

    # 4. Use 'sed' to replace the placeholder with the real absolute path and write to the file.
    #    This is the most robust way to handle paths with special characters.
    echo "\${JSON_TEMPLATE}" | sed "s|___TSV_PATH_PLACEHOLDER___|\${TSV_PATH}|g" > ${meta.id}_wdl_inputs.json
    """
}

process RUN_WDL_VARIANT_CALLING {
    tag "Variant Calling on ${meta.id}"
    publishDir "${params.outdir}/${meta.id}/variant_calling/", mode: 'copy'

    input:
    tuple val(meta), path(wdl_inputs_json)
    path cromwell_config

    output:
    tuple val(meta), path("final_wdl_output"), emit: wdl_results

    tuple val(meta),
          path("final_wdl_output/**/${meta.id}.merged.final.split.vcf"),
          path("final_wdl_output/**/${meta.id}.merged.haplocheck_contamination.txt"),
          path("final_wdl_output/**/${meta.id}.merged.per_base_coverage.tsv"),
          emit: wdl_files

    script:
    """
    #!/bin/bash
    set -e

    # 1. 定义最终输出目录
    FINAL_OUTPUT_DIR="final_wdl_output"
    mkdir -p \${FINAL_OUTPUT_DIR} # 确保目录存在
    
    # 2. 修改 Cromwell Options，让 Cromwell 将输出复制到 FINAL_OUTPUT_DIR
    cat > cromwell_options.json <<EOF
    {
      "final_workflow_outputs_dir": "\${FINAL_OUTPUT_DIR}",
      "final_workflow_log_dir": ".",
      "default_runtime_attributes": {
        "queue": "${params.cromwell_options.queue ?: ''}",
        "cpus": ${params.cromwell_options.cpus},
        "memory": ${params.cromwell_options.memory},
        "runtime_minutes": ${params.cromwell_options.runtime_minutes}
      }
    }
    EOF

    # 3. 运行 Cromwell
    java -Dconfig.file=${cromwell_config} \\
         -jar ${params.cromwell_jar} run \\
         ${params.wdl_script} \\
         --inputs ${wdl_inputs_json} \\
         --options cromwell_options.json
    
    # 4. 确保 Cromwell 成功后，最终输出目录中包含文件
    if [ ! -d "\${FINAL_OUTPUT_DIR}" ] || [ -z "\$(ls -A \${FINAL_OUTPUT_DIR})" ]; then
        echo "ERROR: Cromwell finished but no files were found in the final output directory: \${FINAL_OUTPUT_DIR}" >&2
        exit 1
    fi
    """
}

process ANNOTATE_INDIVIDUAL_VCF {
    tag "Annotate VCF for ${meta.id}"
    publishDir "${params.outdir}/${meta.id}/hail_results/annotation_individual", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(contamination), path(coverage)
    path(hail_dir)

    output:
    tuple val(meta), path("hail_results/annotation_individual/**"), path(".annotate_complete"), emit: annotated_results

    script:
    def config_json = new groovy.json.JsonBuilder(params.hail_pipeline_config).toPrettyString()
    def hail_input_dir_name = "wdl_outputs"

    """
    #!/bin/bash
    set -euo pipefail

    echo "--- [DEBUG] STARTING ANNOTATE_INDIVIDUAL_VCF for ${meta.id} ---"
    echo "[DEBUG] VCF:          ${vcf}"
    echo "[DEBUG] Contamination:${contamination}"
    echo "[DEBUG] Coverage:     ${coverage}"

    # Prepare dirs
    printf '%s' '${config_json}' > config.json
    rm -rf hail_results
    mkdir -p ${hail_input_dir_name}/vcfs 
    mkdir -p ${hail_input_dir_name}/contamination 
    mkdir -p ${hail_input_dir_name}/coverage
    mkdir -p hail_results/annotation_individual/vep_vcf
    mkdir -p hail_results/annotation_individual/metadata
    mkdir -p hail_results/annotation_individual/final_outputs

    # Link files with expected names (no find)
    ln -sf "${vcf}"           "${hail_input_dir_name}/vcfs/${meta.id}.merged.final.split.vcf"
    ln -sf "${contamination}" "${hail_input_dir_name}/contamination/${meta.id}.merged.haplocheck_contamination.txt"
    ln -sf "${coverage}"      "${hail_input_dir_name}/coverage/${meta.id}.merged.per_base_coverage.tsv"

    echo "--- Launching annotation for sample ${meta.id} ---"
    python ${hail_dir}/add_annotation_refine.py --config config.json

    VEP_OUTPUT_FILE="hail_results/annotation_individual/vep_vcf/${meta.id}.merged_vep.vcf"
    if [ -f "\${VEP_OUTPUT_FILE}" ] && [ -s "\${VEP_OUTPUT_FILE}" ]; then
        echo "[SUCCESS] Annotation OK for ${meta.id}"
        touch .annotate_complete
    else
        echo "ERROR: Expected VEP output file '\${VEP_OUTPUT_FILE}' not found or empty." >&2
        echo "--- Debug: listing hail_results ---" >&2
        ls -R hail_results >&2
        exit 1
    fi
    """
}
