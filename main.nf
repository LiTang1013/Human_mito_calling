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

// Create a tuple for the reference genome, including the FASTA and all its indices
ch_ref_genome_tuple = tuple( file(params.ref_genome_fasta), file("${params.ref_genome_fasta}.{amb,ann,bwt,pac,sa,fai}") )
// Create a channel with just the reference FASTA file, needed for CRAM generation
ch_fasta_only = Channel.value( ch_ref_genome_tuple[0] )
// Channel for the Cromwell configuration file
ch_cromwell_conf = file("${baseDir}/cromwell.conf")
// Channel for the directory containing the custom Python/Hail scripts
ch_hail_directory = file("${baseDir}/hail")

// --- Channel for input samples ---
// This channel reads the input TSV (sample_id \t fastq_url),
// groups all URLs by sample_id, sorts them, and then creates
// R1/R2 pairs. Finally, it flattens the list so that each
// FASTQ pair is a separate emission, tagged with a unique 'pair_id'.
// Output: [ [id: 'sample_A', pair_id: 0], 'url_R1', 'url_R2' ]
//         [ [id: 'sample_A', pair_id: 1], 'url_R1_b', 'url_R2_b' ]
ch_samples = Channel.fromPath(params.input)
    .splitCsv(header: false, sep: '\t') 
    .filter { row -> row.size() == 2 && row[0] && !row[0].trim().isEmpty() && row[1] && !row[1].trim().isEmpty() } // Filter malformed rows
    .map { row -> tuple(row[0], row[1]) } // -> [sample_id, url]
    .groupTuple() // -> [sample_id, [url1, url2, ...]]
    .map { sample_id, urls_list ->
        def meta = [id: sample_id]
        def sorted_urls = urls_list.sort() // Sort to ensure R1/R2 alignment (assumes consistent naming)
        if (sorted_urls.size() % 2 != 0) {
            error "Sample ${meta.id} has an odd number of FASTQ files after grouping. Please check input file."
        }
        def fastq_pairs = []
        for (int i = 0; i < sorted_urls.size(); i += 2) {
            fastq_pairs.add([sorted_urls[i], sorted_urls[i+1]]) // Create [R1, R2] pairs
        }
        return [meta, fastq_pairs] // -> [meta, [[R1, R2], [R1_b, R2_b], ...]]
    }
    .flatMap { meta, fastq_pairs ->
        // Emit one item per FASTQ pair, adding a unique pair_id
        fastq_pairs.indexed().collect { i, pair -> [meta + [pair_id: i], pair[0], pair[1]] }
    }

// ====================================================================================
//                                    WORKFLOW
// ====================================================================================
workflow {
    // 1. Download each FASTQ pair
    DOWNLOAD_FASTQ(ch_samples) // Output: [meta (with pair_id), fq1, fq2]

    // 2. Align pair to reference, output unsorted BAM
    ALIGN_AND_UNSORT(DOWNLOAD_FASTQ.out, ch_ref_genome_tuple) // Output: [meta, unsorted_bam]

    // 3. Sort BAM and convert to CRAM for each pair
    SORT_AND_CONVERT_TO_CRAM(ALIGN_AND_UNSORT.out.unsorted_bam, ch_fasta_only) // Output: [meta, single_cram, single_crai]

    // Prepare CRAM files for merging by grouping all pairs belonging to the same sample
    ch_crams_to_merge = SORT_AND_CONVERT_TO_CRAM.out.single_cram
        .map { meta, cram, crai -> [meta.id, [cram, crai]] } // Map to [sample_id, [cram, crai]]
        .groupTuple() // Group by sample_id: [sample_id, [[cram1, crai1], [cram2, crai2], ...]]
        .map { sample_id, cram_pair_list ->
            def meta = [id: sample_id]
            // Collect CRAMs and CRAIs into separate lists
            def crams = cram_pair_list.collect { it[0] }
            def crais = cram_pair_list.collect { it[1] }
            return [meta, crams, crais] // Re-package as [meta, [cram1, ...], [crai1, ...]]
        }

    // 4. Merge all CRAMs for a sample into one final CRAM
    MERGE_CRAMS(ch_crams_to_merge) // Output: [meta, merged_cram, merged_crai]

    // 5. Downstream processes connect to the output of MERGE_CRAMS
    GENERATE_CRAM_TSV(MERGE_CRAMS.out.merged_cram)
    GENERATE_WDL_JSON(GENERATE_CRAM_TSV.out.tsv)
    RUN_WDL_VARIANT_CALLING(GENERATE_WDL_JSON.out.json, ch_cromwell_conf)
    ANNOTATE_INDIVIDUAL_VCF(RUN_WDL_VARIANT_CALLING.out.wdl_files, ch_hail_directory)
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
//                                    PROCESSES
// ====================================================================================

process DOWNLOAD_FASTQ {
    tag "Download for ${meta.id} (Pair ${meta.pair_id})"

    // Publish downloaded FASTQs to results dir. 
    // The local work-dir copy will be cleaned up by ALIGN_AND_UNSORT.
    publishDir "${params.outdir}/${meta.id}/fastq", mode: 'copy' 

    input:
    tuple val(meta), val(url1), val(url2)
    output:
    tuple val(meta), path("${meta.id}_${meta.pair_id}_1.fastq.gz"), path("${meta.id}_${meta.pair_id}_2.fastq.gz")
    script:
    """
    #!/bin/bash
    set -e
    # Wget download
    wget --no-check-certificate -O "${meta.id}_${meta.pair_id}_1.fastq.gz" "${url1}"
    wget --no-check-certificate -O "${meta.id}_${meta.pair_id}_2.fastq.gz" "${url2}"
    """
}


process ALIGN_AND_UNSORT {
    tag "BWA-MEM on ${meta.id} (Pair ${meta.pair_id})"

    // Resource request should match BWA-MEM itself, which is often lower than sorting
    cpus 8 
    memory '32 GB'
    time '1 h'

    input:
    // Receives a single FASTQ pair
    tuple val(meta), path(read1), path(read2)
    tuple path(ref_fasta), path(bwa_indices)
    
    output:
    // Output a single unsorted BAM file
    tuple val(meta), path("${meta.id}_${meta.pair_id}.unsorted.bam"), emit: unsorted_bam 
    
    script:
    // CRITICAL: Read Group ID must be unique. Use meta.id and pair_id combination.
    def read_group_id = "${meta.id}_${meta.pair_id}" 
    // Define the full read group string for BWA
    def read_group = "\'@RG\\tID:${read_group_id}\\tSM:${meta.id}\\tPL:ILLUMINA\'"
    """
    #!/bin/bash
    set -e

    # Run alignment and pipe output to samtools to create an unsorted BAM
    bwa mem -t ${task.cpus} -R ${read_group} ${ref_fasta} ${read1} ${read2} | \\
    samtools view -@ ${task.cpus} -b -o ${meta.id}_${meta.pair_id}.unsorted.bam -
    
    # After successfully creating the BAM, delete the original FASTQ files (local work dir copy only)
    rm ${read1} ${read2}
    """
}

process SORT_AND_CONVERT_TO_CRAM {
    tag "Sort and CRAM for ${meta.id} (Pair ${meta.pair_id})"

    // Increased resources to potentially solve OOM (Out Of Memory) issues during sorting
    cpus 16
    memory '64 GB'
    time '1 h'
    
    input:
    // Receives a single unsorted BAM file
    tuple val(meta), path(unsorted_bam)
    path ref_fasta
    
    output:
    // Output a single CRAM and its index
    tuple val(meta), path("${meta.id}_${meta.pair_id}.cram"), path("${meta.id}_${meta.pair_id}.cram.crai"), emit: single_cram
    
    script:
    // Groovy variable (used only for filename construction, not in Bash logic)
    def unsorted_bam_filename = unsorted_bam.getName()

    """
    #!/bin/bash
    set -e
    
    # 1. Define Bash variables for intermediate and final files
    SORTED_BAM="${meta.id}_${meta.pair_id}.sorted.bam"
    OUTPUT_CRAM="${meta.id}_${meta.pair_id}.cram"

    # 2. Sort: Use samtools sort
    # Note: We use the Nextflow-provided filename ${unsorted_bam_filename} here
    samtools sort -@ ${task.cpus} -o \${SORTED_BAM} ${unsorted_bam_filename}
    
    # 3. Convert to CRAM format, using the reference fasta
    samtools view -@ ${task.cpus} -T ${ref_fasta} -C -o \${OUTPUT_CRAM} \${SORTED_BAM}

    # 4. Create CRAM index
    samtools index \${OUTPUT_CRAM}

    # 5. Cleanup: Delete intermediate files (unsorted and sorted BAMs)
    rm ${unsorted_bam_filename} \${SORTED_BAM}
    """
}

process MERGE_CRAMS {
    tag "Merge CRAMs for ${meta.id}"

    // Publish the final merged CRAM file for the sample
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy', pattern: "*.{cram,crai}"

    input:
    // Receives a list of all CRAMs/CRAIs for this sample
    tuple val(meta), path(crams), path(crais)

    output:
    // Output the merged CRAM and CRAI
    tuple val(meta), path("${meta.id}.merged.cram"), path("${meta.id}.merged.cram.crai"), emit: merged_cram

    script:
    // Join the file lists into space-separated strings for the command line
    def cram_files = crams.join(' ')
    def crai_files = crais.join(' ')

    """
    #!/bin/bash
    set -e
    
    # 1. Define Bash variables, ensure output filename is correct
    OUTPUT_CRAM="${meta.id}.merged.cram"

    # 2. Use samtools merge to combine all CRAM files
    # Note: -f forces overwrite, -O CRAM specifies output format
    samtools merge -@ ${task.cpus} -f -O CRAM --output-fmt-option version=3.0 -o \${OUTPUT_CRAM} ${cram_files}

    # 3. Create CRAM index for the merged file
    samtools index \${OUTPUT_CRAM}
    
    # 4. Cleanup: Delete all original single-pair CRAM files and their indices
    rm ${cram_files} ${crai_files}
    """
}

// ... (GENERATE_CRAM_TSV, GENERATE_WDL_JSON, RUN_WDL_VARIANT_CALLING, ANNOTATE_INDIVIDUAL_VCF processes remain unchanged)

process GENERATE_CRAM_TSV {
    tag "Generate CRAM TSV for ${meta.id}"

    // This TSV file is an input for the WDL pipeline
    publishDir "${params.outdir}/${meta.id}/variant_calling/inputs", mode: 'copy'
    input:
    tuple val(meta), path(cram), path(crai)
    output:
    tuple val(meta), path("${meta.id}_cram_list.tsv"), emit: tsv
    script:
    """
    # Get the full, absolute path of the files, required by WDL
    cram_path=\$(readlink -f ${cram})
    crai_path=\$(readlink -f ${crai})
    # Create a tab-separated file listing the CRAM and its index
    echo -e "\${cram_path}\\t\${crai_path}" > ${meta.id}_cram_list.tsv
    """
}


process GENERATE_WDL_JSON {
    tag "Generate WDL JSON for ${meta.id}"
    // This JSON file is the main input config for Cromwell
    publishDir "${params.outdir}/${meta.id}/variant_calling/inputs", mode: 'copy'

    input:
    tuple val(meta), path(cram_tsv)

    output:
    tuple val(meta), path("${meta.id}_wdl_inputs.json"), emit: json

    script:
    // 1. In Groovy, create a JSON template with a unique placeholder.
    def json_content = new groovy.json.JsonBuilder()
    // Add the WDL-specific pipeline namespace to each parameter key
    def wdl_inputs = params.wdl_inputs.collectEntries { key, value ->
        ["MitochondriaMultiSamplePipeline.${key}", value]
    }
    // The inputSamplesFile must be the absolute path generated in the bash script
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

    // We define two 'output' blocks to capture different sets of files
    // 'wdl_results' captures the entire output directory (for debugging/archive)
    output:
    tuple val(meta), path("final_wdl_output"), emit: wdl_results

    // 'wdl_files' captures the specific files needed for the next annotation step
    output:
    tuple val(meta),
      path("final_wdl_output/**/execution/${meta.id}.merged.final.split.vcf"),
      path("final_wdl_output/**/execution/${meta.id}.merged.haplocheck_contamination.txt"),
      path("final_wdl_output/**/execution/${meta.id}.merged.per_base_coverage.tsv"),
      emit: wdl_files

    script:
    """
    #!/bin/bash
    set -e

    # 1. Define the final output directory
    FINAL_OUTPUT_DIR="final_wdl_output"
    mkdir -p \${FINAL_OUTPUT_DIR} # Ensure directory exists
    
    # 2. Create Cromwell Options JSON to copy final outputs to FINAL_OUTPUT_DIR
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

    # 3. Run Cromwell
    java -Dconfig.file=${cromwell_config} \\
         -jar ${params.cromwell_jar} run \\
         ${params.wdl_script} \\
         --inputs ${wdl_inputs_json} \\
         --options cromwell_options.json
    
    # 4. Verify that Cromwell succeeded and files exist in the output directory
    if [ ! -d "\${FINAL_OUTPUT_DIR}" ] || [ -z "\$(ls -A \${FINAL_OUTPUT_DIR})" ]; then
        echo "ERROR: Cromwell finished but no files were found in the final output directory: \${FINAL_OUTPUT_DIR}" >&2
        exit 1
    fi
    """
}

process ANNOTATE_INDIVIDUAL_VCF {
    tag "Annotate VCF for ${meta.id}_rerun3"
    publishDir "${params.outdir}/${meta.id}/hail_results/annotation_individual", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(contamination), path(coverage)
    path(hail_dir) // This is the directory containing 'add_annotation_refine.py'

    output:
    // Emit all results, plus a '.annotate_complete' file as a success flag
    tuple val(meta), path("hail_results/annotation_individual/**"), path(".annotate_complete"), emit: annotated_results

    script:
    // Create the config.json for the python script from 'params'
    def config_json = new groovy.json.JsonBuilder(params.hail_pipeline_config).toPrettyString()
    // Define the input directory name that the python script expects
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
    # Clean up any previous run and create the expected directory structure
    rm -rf hail_results
    mkdir -p ${hail_input_dir_name}/vcfs 
    mkdir -p ${hail_input_dir_name}/contamination 
    mkdir -p ${hail_input_dir_name}/coverage
    mkdir -p hail_results/annotation_individual/vep_vcf
    mkdir -p hail_results/annotation_individual/metadata
    mkdir -p hail_results/annotation_individual/final_outputs

    # Copy files into the expected directory structure with the expected names
    cp "${vcf}"           "${hail_input_dir_name}/vcfs/${meta.id}.merged.final.split.vcf"
    cp "${contamination}" "${hail_input_dir_name}/contamination/${meta.id}.merged.haplocheck_contamination.txt"
    cp "${coverage}"      "${hail_input_dir_name}/coverage/${meta.id}.merged.per_base_coverage.tsv"

    echo "--- Launching annotation for sample ${meta.id} ---"
    # Launch the main python annotation script
    python ${hail_dir}/add_annotation_refine.py --config config.json

    # Define the expected final output file for validation
    FINAL_OUTPUT_FILE="hail_results/annotation_individual/final_outputs/POC_variant_list.txt"

    # Check if the final file exists (-f) and is not empty (-s)
    if [ -f "\${FINAL_OUTPUT_FILE}" ] && [ -s "\${FINAL_OUTPUT_FILE}" ]; then
        echo "[SUCCESS] Annotation OK for ${meta.id}. Final output '\${FINAL_OUTPUT_FILE}' found."
        # If the file is valid, touch the success flag
        touch .annotate_complete
    else
        # If file is missing or empty, print error and exit with an error code
        echo "ERROR: Expected FINAL output file '\${FINAL_OUTPUT_FILE}' not found or empty." >&2
        echo "This likely means the python script failed silently." >&2
        echo "--- Debug: listing hail_results tree ---" >&2
        ls -R hail_results >&2
        # This ensures Nextflow reports the task as failed
        exit 1
    fi
    """
}