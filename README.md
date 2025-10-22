# Human mtDNA Variant Calling Pipeline
This Nextflow pipeline is designed to perform human mitochondrial (mtDNA) variant calling, starting from raw FASTQ files or URLs. It utilizes the Broad Institute's robust mt-DNA variant calling workflow (implemented in WDL and executed via Cromwell) and a subsequent Hail-based annotation workflow.
The pipeline is primarily configured to run on a Slurm cluster environment using a three-stage approach managed by an SLURM Job Array Manager.

## 1. Prerequisites
Before running the pipeline, ensure you have the following installed and configured:
* Nextflow (v24.10.2 or later recommended)
* Conda/Miniconda (You will need a dedicated conda environment. We recommend using our install script to create the `bio_pipeline_env`)
* Slurm Workload Manager (required for the cluster profile)
* Java Runtime Environment (JRE) (required for Cromwell, Picard, and GATK)
* VEP (v113.3 or later) (required for annotation)

All required software JARs and reference files must be accessible at the absolute paths specified in nextflow.config.

## 2. Pipeline Overview
### Stage 1: Preprocessing (Per-Sample, Per-Pair)

This stage prepares the raw data for variant calling and is executed as part of the Slurm Job Array.

 `DOWNLOAD_FASTQ:` Downloads FASTQ files from provided URLs.

 `ALIGN_AND_UNSORT:` Aligns FASTQ pairs to the reference genome using BWA-MEM and outputs an unsorted BAM file.

 `SORT_AND_CONVERT_TO_CRAM:` Sorts the BAM, converts it to CRAM format (using the reference FASTA), and indexes it.

 `MERGE_CRAMS:` Collects and merges all CRAM files belonging to the same sample into a single merged CRAM file and index.

### Stage 2: Variant Calling and Initial Annotation (Per-Sample)

This stage executes the core WDL workflow and performs initial per-sample annotation, running sequentially after the merging in the Job Array.

 `GENERATE_CRAM_TSV:` Creates a TSV file listing the final merged CRAM and its index.

 `GENERATE_WDL_JSON:` Creates the WDL input JSON, incorporating the sample-specific CRAM TSV path and all static reference/software paths from nextflow.config.

 `RUN_WDL_VARIANT_CALLING:` Executes the external WDL script (MitochondriaPipeline_multisample.wdl) via Cromwell in Slurm mode to perform the core variant calling.
 
 `ANNOTATE_INDIVIDUAL_VCF:` Runs the initial Hail-based annotation on the single-sample VCF output by the WDL workflow.

### Stage 3: Final Merging and Annotation (All Samples)

This stage is executed as a single, dependent Slurm job (--finalize mode) after all Stage 1 and Stage 2 jobs successfully complete. It uses custom Python/Hail scripts to:

* Create a unified Hail Matrix Table for per-base coverage across all samples
* Merge all single-sample VCFs (*.final.vcf.gz) into a single, combined Matrix Table.
* Add final, comprehensive annotations using multiple public and internal databases defined in params.hail_pipeline_config.


## 3. Configuration

All essential parameters and resource settings are managed across two files.

### 3.1. nextflow.config (Pipeline Parameters and Profiles)
This file controls Nextflow's behavior, I/O paths, software paths, and resource defaults.
<div style="overflow-x: auto;">

| Parameter Group   | Key                                      | Description                                                                                                                                                                                                                     | Example Value                |
|--------------------|-------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------|
| **I/O**            | `params.input`                           | Path to the **Master Sample List TSV** file (**Required**).                                                                                                                                                                     | `/path/to/test_sample.tsv`   |
|                    | `params.outdir`                          | Root directory for all pipeline results.                                                                                                                                                                                        | `/path/to/results`           |
|                    | `params.ref_genome_fasta`                | Path to the full human reference FASTA (e.g., *Homo_sapiens_assembly38.fasta*).                                                                                                                                                 | `/path/to/ref/h38.fasta`     |
| **Core Software**  | `params.wdl_script`                      | Path to the WDL script file (e.g., *MitochondriaPipeline_multisample.wdl*).                                                                                                                                                     | `/path/to/script.wdl`        |
|                    | `params.cromwell_jar`                    | Path to the Cromwell JAR file.                                                                                                                                                                                                  | `/path/to/cromwell-90.jar`   |
|                    |                                           | A map containing all static paths needed by the WDL (Picard, GATK, specific reference files like chrM fasta/dict/indices, blacklist BEDs, etc.). **These must all be absolute paths.**                                          | *(See file content)*         |
| **WDL Inputs**     | `params.wdl_inputs`                      | A map for all WDL input settings.                                                                                                                                                                                               | *(See file content)*         |
| **Hail Config**    | `params.hail_pipeline_config`             | A map for Hail-related settings, including VEP cache directory, input metadata, and all annotation database paths. **These must all be absolute paths.**                                                                        | *(See file content)*         |
| **Profiles**       | `profiles.cluster`                       | Uses the **slurm** executor and enables Conda (environment: `/home/lt692/.conda/envs/bio_pipeline_env`). It sets the internal WDL engine's backend to **Slurm** via the `cromwell.backend` setting.                            |                              |
| **Process Resources** | `process.withName: '...'`              | Specific resource requests (`cpus`, `memory`, `time`) for heavy processes like **ALIGN_AND_SORT** (16 cpus, 18 GB) and **SORT_AND_CONVERT_TO_CRAM** (16 cpus, 64 GB).                                                           | *(See file content)*         |

</div>


### 3.2. launch_pipeline.sh (Slurm Manager Script)
This is the primary script used to launch the three-stage process on the Slurm cluster.

<div style="overflow-x: auto;">

| Variable | Description | Default Path (Needs Review) |
|-----------|--------------|------------------------------|
| **MAIN_NF_SCRIPT** | Path to the Nextflow pipeline script (*main.nf*). | `/path/to/main.nf` |
| **NEXTFLOW_CONFIG** | Path to the configuration file (*nextflow.config*). | `/path/to/nextflow.config` |
| **MASTER_SAMPLE_LIST** | **The only path you should modify here;** this must match `params.input` in *nextflow.config*. | `/path/to/test_sample.tsv` |
| **HAIL_SCRIPT_DIR** | Directory containing the Hail/Python scripts for Stage 2. | `/path/to/hail` |
| **OUTPUT_DIR** | Root output directory (must match `params.outdir` in *nextflow.config*). | `/path/to/results` |
| **CONCURRENT_SAMPLES** | The maximum number of array tasks (samples) to run simultaneously in Stage 1. | `3` |

</div>

## 4. Input Format
The pipeline requires a Master Sample List TSV file (e.g., test_sample.tsv).

Format: A tab-separated file with two columns, no header.
```
Sample1	https://example.com/data/Sample1_L1_R1.fastq.gz
Sample1	https://example.com/data/Sample1_L1_R2.fastq.gz
Sample2	https://example.com/data/Sample2_L1_R1.fastq.gz
Sample2	https://example.com/data/Sample2_L1_R2.fastq.gz
Sample2	https://example.com/data/Sample2_L2_R1.fastq.gz
Sample2	https://example.com/data/Sample2_L2_R2.fastq.gz
```

## 5. Execution
The entire pipeline is launched using the provided Slurm submission script. Do not run nextflow run directly unless you are debugging a single step.

Review and Update: Ensure all paths in launch_pipeline.sh and nextflow.config are correct for your system.

Launch the Pipeline:
```bash
bash launch_pipeline.sh
```

The launch_pipeline.sh script will perform the following automatically:

Master Mode: Calculate the total number of samples and submit the Stage 1 Job Array.

Worker Mode (Stage 1): Each array task will run the Nextflow pipeline in a dedicated work directory (-w) for one sample, using the cluster profile.

Finalizer Mode (Stage 2): A final job is submitted with an afterok dependency on the array, executing the Hail merging and final annotation only after all sample-specific jobs succeed.

## 6. Output Structure
All results are deposited in the directory specified by params.outdir (/path/to/results in the example).
The structure is:
```
/path/to/results/
├── Sample1/
│   ├── fastq/
│   │   ├── Sample1_0_1.fastq.gz
│   │   └── Sample1_0_2.fastq.gz
│   ├── alignment/
│   │   ├── Sample1.merged.cram
│   │   └── Sample1.merged.cram.crai
│   ├── variant_calling/
│   │   ├── inputs/
│   │   │   ├── Sample1_cram_list.tsv
│   │   │   └── Sample1_wdl_inputs.json
│   │   └── final_wdl_output/ # Output from WDL/Cromwell
│   │       ├── Sample1.final.vcf.gz
│   │       ├── Sample1.per_base_coverage.tsv
│   │       └── ... (other WDL outputs)
│   └── hail_results/
│       └── annotation_individual/ # Stage 1 Hail output
│           ├── Sample1.final.annotated.vcf.gz
│           └── ...
└── (Root Level - Stage 2 Hail Final Outputs)
    ├── coverage.mt                 # Hail Matrix Table for coverage
    ├── combined_variants.mt        # Merged, filtered VCF Matrix Table
    └── final_annotated_results.mt  # Final annotated Matrix Table (primary output)
```
