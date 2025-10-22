import os
import sys
import subprocess
import glob
import pandas as pd
import io
import json
import csv
from concurrent.futures import ThreadPoolExecutor
import functools

# ==============================================================================
# Helper Functions
# ==============================================================================

def run_command(command, cwd=None, is_shell=False):
    """
    Executes a command and streams its output in real-time.
    """
    cmd_str = command[0] if is_shell else ' '.join(command)
    print(f"[*] Executing (in {cwd or '.'}): {cmd_str}")
    
    # Use Popen to stream output in real-time
    process = subprocess.Popen(
        cmd_str,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=cwd,
        shell=is_shell,
        executable='/bin/bash'
    )
    
    # Read and print output line by line
    for line in process.stdout:
        print(line, end='')
    
    # Wait for the process to complete and check the return code
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd_str)
    
    print(f"[+] Command successful: {cmd_str}")


def read_files_to_dataframe(pattern, func_extract, col_names):
    all_data = []
    for filepath in glob.glob(pattern):
        with open(filepath, 'r') as f:
            df = pd.read_csv(f, sep='\t')
            all_data.append(func_extract(df))
    if not all_data: return pd.DataFrame(columns=col_names)
    combined_df = pd.concat(all_data, ignore_index=True)
    combined_df.columns = col_names
    return combined_df

def _run_vep_task(file_path, config):
    """
    Helper function to run VEP for a single file. Designed for parallel execution.
    """
    try:
        output_dir = config['vep_vcf_dir']
        cache_dir = config['vep_cache_dir']
        
        filename = os.path.basename(file_path).replace('.final.split.vcf', '')
        output_path = os.path.join(output_dir, f"{filename}_vep.vcf")
        
        print(f"[*] Starting VEP for: {filename}")
        
        vep_args = ['vep', '--cache', '-i', file_path, '-o', output_path, '--dir_cache', cache_dir, '--no_stats', '--distance', '0', '--biotype', '--symbol', '--hgvs', '--variant_class', '--force_overwrite', '--vcf']
        vep_command_str = ' '.join(vep_args)
        full_command_string = f"module load VEP/113.3-GCC-13.3.0 && {vep_command_str}"
        
        run_command([full_command_string], is_shell=True)
        print(f"[+] Finished VEP for: {filename}")
        return f"Success: {filename}"
    except Exception as e:
        print(f"[!] ERROR processing {file_path}: {e}")
        return f"Error: {os.path.basename(file_path)}"
        
# --- MODIFICATION START: Added helper for Picard ---
def _run_picard_task(cram_path, config):
    """Helper function to run Picard CollectWgsMetrics for a single CRAM file."""
    try:
        output_dir = config['wgs_metrics_dir']
        ref_genome = config['reference_genome_path']
        
        sample_name = os.path.basename(cram_path).split('.')[0]
        output_path = os.path.join(output_dir, f"{sample_name}.wgsMetrics.txt")
        
        print(f"[*] Starting Picard CollectWgsMetrics for: {sample_name}")
        
        # This command assumes picard is in the path after 'module load'.
        # $EBROOTPICARD is a common environment variable on HPCs for the jar path.
        picard_command = (
            "java -Xmx4g -jar $EBROOTPICARD/picard.jar CollectWgsMetrics "
            f"I={cram_path} "
            f"O={output_path} "
            f"R={ref_genome}"
        )
        
        full_command_string = f"module load picard && {picard_command}"
        
        run_command([full_command_string], is_shell=True)
        print(f"[+] Finished Picard for: {sample_name}")
        return f"Success: {sample_name}"
    except Exception as e:
        print(f"[!] ERROR processing {cram_path}: {e}")
        return f"Error: {os.path.basename(cram_path)}"
# --- MODIFICATION END ---


# ==============================================================================
# Main Analysis Steps
# ==============================================================================
def step1_run_vep(config):
    print("\n--- Step 1: Running VEP Annotation in Parallel ---")
    
    input_dir, output_dir = config['raw_vcf_dir'], config['vep_vcf_dir']
    os.makedirs(output_dir, exist_ok=True)
    
    vcf_files = glob.glob(os.path.join(input_dir, '*.final.split.vcf'))
    if not vcf_files:
        print(f"[!] No VCF files found in {input_dir}, skipping VEP.")
        return

    num_workers = config.get('num_workers', os.cpu_count())
    print(f"[*] Starting VEP annotation for {len(vcf_files)} files using up to {num_workers} parallel threads...")

    task_with_config = functools.partial(_run_vep_task, config=config)

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(task_with_config, vcf_files))

    success_count = sum(1 for r in results if r.startswith("Success"))
    error_count = len(results) - success_count
    print(f"\n[+] VEP annotation complete. {success_count} files processed successfully, {error_count} failed.")

def step2_extract_metadata(config):
    print("\n--- Step 2: Extracting Haplogroup and Contamination Info ---")
    source_dir, output_dir = config['haplocheck_dir'], config['metadata_dir']
    os.makedirs(output_dir, exist_ok=True)
    file_pattern = os.path.join(source_dir, "*.txt")
    haplogroup_df = read_files_to_dataframe(file_pattern, lambda df: df[['SampleID', 'HgMajor']], ['Sample_ID', 'Full_Haplogroup'])
    hap_path = os.path.join(output_dir, "haplogroup_full.txt")
    haplogroup_df.to_csv(hap_path, sep='\t', index=False)
    print(f"[+] Haplogroup file saved to: {hap_path}")
    contam_df = read_files_to_dataframe(file_pattern, lambda df: df[['SampleID', 'Contamination']], ['Sample_ID', 'Contamination_Status'])
    contam_path = os.path.join(output_dir, "contamination.txt")
    contam_df.to_csv(contam_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)
    print(f"[+] Contamination file saved to: {contam_path}")
    print(f"[+] Sample contamination statistics:\n{contam_df['Contamination_Status'].value_counts().to_string()}")

def step3_run_variant_processor(config):
    """Calls the external script to process and filter variants."""
    print("\n--- Step 3: Calling Variant Processing Script ---")
    
    python_command_list = [
        sys.executable, "process_variants.py",
        "--vep-vcf-dir", config["vep_vcf_dir"],
        "--final-output-dir", config["final_output_dir"],
        "--fullhaplogroups", config["fullhaplogroups"],
        "--contamination", config["contamination"],
        "--master-pedigree-file", config["master_pedigree_file"],
        "--gnomadcache", config["gnomadcache"],
        "--clinvarcache", config["clinvarcache"],
        "--mitomap-polycache", config["mitomap_polycache"],
        "--mitomap-diseasecache", config["mitomap_diseasecache"],
        "--helixcache", config["helixcache"],
        "--haplogroup-varcache", config["haplogroup_varcache"],
        "--mitimpactcache", config["mitimpactcache"],
        "--mitotipcache", config["mitotipcache"],
        "--hmtvarcache", config["hmtvarcache"],
    ]
    python_command_str = ' '.join(python_command_list)
    run_command([python_command_str], is_shell=True)

# --- MODIFICATION START: Added new step 4 ---
def step4_generate_wgs_metrics(config):
    """
    NEW STEP: Generate WGS metrics files from CRAM files using Picard.
    """
    print("\n--- Step 4: Generating WGS Metrics from CRAM files ---")
    output_dir = config['wgs_metrics_dir']
    os.makedirs(output_dir, exist_ok=True)

    try:
        # Assumes a single column TSV with full paths to CRAM files
        samples_df = pd.read_csv(config['samples_tsv_path'], sep='\t', header=None)
        cram_files = samples_df.iloc[:, 0].tolist()
    except FileNotFoundError:
        print(f"[!] ERROR: Samples TSV file not found at {config['samples_tsv_path']}. Cannot generate WGS metrics.")
        return
    
    if not cram_files:
        print(f"[!] No CRAM files found in {config['samples_tsv_path']}, skipping WGS metrics generation.")
        return

    num_workers = config.get('num_workers', os.cpu_count())
    print(f"[*] Starting Picard CollectWgsMetrics for {len(cram_files)} files using up to {num_workers} parallel threads...")

    task_with_config = functools.partial(_run_picard_task, config=config)

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(task_with_config, cram_files))

    success_count = sum(1 for r in results if r.startswith("Success"))
    error_count = len(results) - success_count
    print(f"\n[+] Picard WGS metrics generation complete. {success_count} files processed successfully, {error_count} failed.")
# --- MODIFICATION END ---

def step5_calculate_nuclear_coverage(config):
    print("\n--- Step 5 (was 4): Calculating Median Nuclear Coverage ---")
    metrics_files = glob.glob(os.path.join(config['wgs_metrics_dir'], '*wgsMetrics.txt'))
    print(f"--- [DIAGNOSTIC] Found {len(metrics_files)} files matching '*wgsMetrics.txt' in {config['wgs_metrics_dir']}")
    
    coverage_data = []

    for file_path in metrics_files:
        print(f"--- [DIAGNOSTIC] Processing file: {os.path.basename(file_path)}")
        sample = os.path.basename(file_path).replace('.wgsMetrics.txt', '')
        # --- MODIFICATION START: Correctly parse Picard output ---
        # The data is under a header and separated by a blank line
        with open(file_path, 'r') as f:
            lines = f.readlines()
            header_index = -1
            for i, line in enumerate(lines):
                if line.startswith('MEDIAN_COVERAGE'):
                    header_index = i
                    break
            
            if header_index != -1 and header_index + 2 < len(lines):
                data_line = lines[header_index + 2]
                median_cov = data_line.strip().split('\t')[1]
                coverage_data.append({'sample': sample, 'median_nuc_cov': median_cov})
    # --- MODIFICATION END ---
    
    if not coverage_data:
        print("--- [DIAGNOSTIC] WARNING: No median coverage data was extracted. The output file will be empty.")
        
    nuc_cov_df = pd.DataFrame(coverage_data)
    output_path = os.path.join(config['final_output_dir'], f"nuc_median_cov_all.txt")
    nuc_cov_df.to_csv(output_path, sep='\t', index=False)
    print(f"[+] Median nuclear coverage file saved to: {output_path}")
    return output_path

def step6_calculate_mtcn(config, nuc_cov_path):
    print("\n--- Step 6 (was 5): Calculating mtDNA Copy Number ---")
    cov_files = glob.glob(os.path.join(config['mt_coverage_dir'], '*.tsv'))
    mean_cov_data = []
    for file_path in cov_files:
        sample_id = os.path.basename(file_path).split('.')[0]
        df = pd.read_csv(file_path, sep='\t'); mean_cov = df['coverage'].mean()
        mean_cov_data.append({'sample': sample_id, 'mean_coverage': mean_cov})
    mean_coverage_df = pd.DataFrame(mean_cov_data)

    if not os.path.exists(nuc_cov_path) or os.path.getsize(nuc_cov_path) == 0:
        print(f"[!] WARNING: Nuclear coverage file is missing or empty at {nuc_cov_path}. Cannot calculate mtCN.")
        # Create an empty file to avoid crashing later steps and provide a clear output
        output_path = os.path.join(config['final_output_dir'], f"mtDNA_CN.txt")
        pd.DataFrame(columns=['Sample_ID', 'Mean_mtDNA_Coverage', 'Median_Nuclear_Coverage', 'mtCopyNumber']).to_csv(output_path, sep='\t', index=False)
        print(f"[+] Empty mtDNA copy number file created at: {output_path}")
        return

    median_df = pd.read_csv(nuc_cov_path, sep='\t', dtype=str)
    # The sample name from Picard might have suffixes, so we match on the start of the string
    mean_coverage_df['join_key'] = mean_coverage_df['sample'].str.split('-').str[0]
    median_df['join_key'] = median_df['sample'].str.split('-').str[0]

    merged_data = pd.merge(mean_coverage_df, median_df, on='join_key', how='left')
    merged_data.rename(columns={'sample_x': 'Sample_ID', 'mean_coverage': 'Mean_mtDNA_Coverage', 'median_nuc_cov': 'Median_Nuclear_Coverage'}, inplace=True)
    merged_data['Mean_mtDNA_Coverage'] = pd.to_numeric(merged_data['Mean_mtDNA_Coverage'], errors='coerce')
    merged_data['Median_Nuclear_Coverage'] = pd.to_numeric(merged_data['Median_Nuclear_Coverage'], errors='coerce')
    merged_data['mtCopyNumber'] = 2 * (merged_data['Mean_mtDNA_Coverage'] / merged_data['Median_Nuclear_Coverage'])
    output_path = os.path.join(config['final_output_dir'], f"mtDNA_CN.txt")
    
    # Select and reorder final columns
    final_cols = ['Sample_ID', 'Mean_mtDNA_Coverage', 'Median_Nuclear_Coverage', 'mtCopyNumber']
    merged_data[final_cols].to_csv(output_path, sep='\t', index=False, na_rep='NA')
    print(f"[+] mtDNA copy number file saved to: {output_path}")


# ==============================================================================
# Main Execution Block
# ==============================================================================
def main():
    with open('config.json', 'r') as f:
        config = json.load(f)
    
    os.makedirs(config['vep_vcf_dir'], exist_ok=True)
    os.makedirs(config['metadata_dir'], exist_ok=True)
    os.makedirs(config['final_output_dir'], exist_ok=True)

    try:
        step1_run_vep(config)
        step2_extract_metadata(config)
        step3_run_variant_processor(config)
        # --- MODIFICATION START: Call the new steps in order ---
        step4_generate_wgs_metrics(config)
        nuc_cov_path = step5_calculate_nuclear_coverage(config)
        step6_calculate_mtcn(config, nuc_cov_path)
        # --- MODIFICATION END ---
        
        print("\n[SUCCESS] All steps in the pipeline completed successfully!")

    except Exception as e:
        print(f"\n[ERROR] The pipeline was terminated due to an error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

