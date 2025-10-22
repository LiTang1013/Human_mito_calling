#!/usr/bin/env python
import os
import glob
import pandas as pd
import io
import json
import argparse

# ==============================================================================
# Helper Function
# ==============================================================================

def load_complex_db(path, key_cols, val_cols, delimiter='\t'):
    """A robust helper function to load various annotation files into dictionaries."""
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment='#')
    index = pd.MultiIndex.from_frame(df[key_cols])
    df['__values__'] = df[val_cols].apply(tuple, axis=1)
    return pd.Series(df['__values__'].values, index=index).to_dict()

# ==============================================================================
# Main Processing Logic
# ==============================================================================

def process_and_filter_variants(args):
    """
    Main logic for processing, annotating, and filtering variants.
    This function contains the logic from the original step3.
    """
    print("\n--- Starting Variant Processing, Annotation, and Filtering ---")
    
    print("[*] Loading metadata and all annotation databases...")
    hap_dict = pd.read_csv(args.fullhaplogroups, sep='\t', header=0, dtype=str, index_col=0).squeeze().to_dict()
    contam_dict = pd.read_csv(args.contamination, sep='\t', header=0, dtype=str, index_col=0).squeeze().to_dict()
    master_pedigree_df = pd.read_csv(args.master_pedigree_file, sep='\t', header=0, dtype=str)
    category_dict = pd.Series(master_pedigree_df.category.values, index=master_pedigree_df.sample_id).to_dict()
    
    # Load standard annotation databases
    gnomad_db = load_complex_db(args.gnomadcache, ['ref', 'position', 'alt'], ['max_observed_heteroplasmy', 'AF_hom', 'AF_het'])
    mitomap_poly_db = load_complex_db(args.mitomap_polycache, ['ref', 'pos', 'alt'], ['gbcnt'])
    mitomap_disease_db = load_complex_db(args.mitomap_diseasecache, ['ref', 'pos', 'alt'], ['status', 'homoplasmy', 'heteroplasmy', 'disease'])
    haplovar_df = pd.read_csv(args.haplogroup_varcache, sep='\t', dtype=str)
    haplovar_db = pd.Series(haplovar_df.Assoc_haplogroups.values, index=haplovar_df.Variant).str.lower().to_dict()
    apogee_db = load_complex_db(args.mitimpactcache, ['Ref', 'Start', 'Alt'], ['APOGEE1', 'APOGEE2'])
    hmtvar_db = load_complex_db(args.hmtvarcache, ['REF', 'POS', 'ALT'], ['HmtVar'])
    mitotip_map = {'Q1': "likely pathogenic", 'Q2': "possibly pathogenic", 'Q3': "possibly benign", 'Q4': "likely benign"}
    mitotip_df = pd.read_csv(args.mitotipcache, sep='\t', dtype=str)
    mitotip_df['prediction'] = mitotip_df['Quartile'].map(mitotip_map)
    mitotip_db = pd.Series(mitotip_df.prediction.values, index=tuple(zip(mitotip_df['rCRS'], mitotip_df['Position'], mitotip_df['Alt']))).to_dict()

    # Special handling for HelixMTdb
    print("[*] Parsing complex HelixMTdb format...")
    helix_df = pd.read_csv(args.helixcache, sep='\t', dtype=str)
    helix_df = helix_df[helix_df["alleles"].str.count(",") == 1].copy()
    helix_df['pos'] = helix_df["locus"].str.split("chrM:").str[1]
    alleles_split = helix_df["alleles"].str.split('"').str
    helix_df['ref'] = alleles_split[1]
    helix_df['alt'] = alleles_split[3]
    helix_df[['AF_hom', 'max_ARF']] = helix_df[['AF_hom', 'max_ARF']].astype(float)
    helix_df['max_het'] = helix_df.apply(lambda row: 1.0 if row['AF_hom'] > 0 else row['max_ARF'], axis=1)
    helix_db_index = pd.MultiIndex.from_frame(helix_df[['ref', 'pos', 'alt']])
    helix_db_series = pd.Series(zip(helix_df['max_het'], helix_df['AF_hom'], helix_df['AF_het']), index=helix_db_index)
    helix_db = helix_db_series.to_dict()
    
    # Special handling for ClinVar
    print("[*] Parsing complex ClinVar format...")
    clinvar_df = pd.read_csv(args.clinvarcache, sep='\t', dtype=str)
    clinvar_df['pos'] = clinvar_df['GRCh38Location']
    spdi_split = clinvar_df['Canonical SPDI'].str.split(':', expand=True)
    clinvar_df['ref'] = spdi_split[2]
    clinvar_df['alt'] = spdi_split[3]
    clinvar_df = clinvar_df[(clinvar_df['ref'].str.len() == 1) & (clinvar_df['alt'].str.len() == 1) & (clinvar_df['ref'] != clinvar_df['alt']) & (clinvar_df['Germline classification'] != "Conflicting interpretations of pathogenicity")].copy()
    clinvar_db_index = pd.MultiIndex.from_frame(clinvar_df[['ref', 'pos', 'alt']])
    clinvar_db_series = pd.Series(clinvar_df['Germline classification'].values, index=clinvar_db_index)
    clinvar_db = clinvar_db_series.to_dict()

    print("[*] Processing and annotating VEP VCF files...")
    all_variants_list = []
    vcf_files = glob.glob(os.path.join(args.vep_vcf_dir, '*_vep.vcf'))
    
    for file_path in vcf_files:
        with open(file_path, 'r') as f: lines = [l for l in f if not l.startswith('##')]
        vcf_df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', dtype=str).rename(columns={'#CHROM': 'CHROM'})
        if vcf_df.empty: continue
        sample_id_col_name = vcf_df.columns[9]
        vcf_df.rename(columns={sample_id_col_name: 'SAMPLE_DATA'}, inplace=True)
        
        if contam_dict.get(sample_id_col_name) == "YES": continue
        vcf_df = vcf_df[~vcf_df['FILTER'].str.contains('base_qual|strand_bias|weak_evidence|blacklisted_site|contamination|position')]
        if vcf_df.empty: continue

        vcf_df[['GT', 'AD', 'AF', 'DP', 'Other']] = vcf_df['SAMPLE_DATA'].str.split(':', n=4, expand=True)
        info_cols = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'HGVS_OFFSET']
        info_split = vcf_df['INFO'].str.split('|', expand=True)
        vcf_df[info_cols] = info_split.iloc[:, 1:len(info_cols)+1]

        vcf_df['SAMPLE_ID'] = sample_id_col_name
        vcf_df['Haplogroup'] = vcf_df['SAMPLE_ID'].map(hap_dict)
        vcf_df['Sample_Category'] = vcf_df['SAMPLE_ID'].str.split('-').str[0].map(category_dict)
        
        var_tuples = list(zip(vcf_df['REF'], vcf_df['POS'], vcf_df['ALT']))
        vcf_df['gnomad_max_hl'], vcf_df['gnomad_af_hom'], vcf_df['gnomad_af_het'] = zip(*[gnomad_db.get(vt, ('0', '0', '0')) for vt in var_tuples])
        vcf_df['helix_max_hl'], vcf_df['helix_af_hom'], vcf_df['helix_af_het'] = zip(*[helix_db.get(vt, (0.0, 0.0, 0.0)) for vt in var_tuples])
        vcf_df['clinvar_interp'] = [clinvar_db.get(vt, '') for vt in var_tuples]
        vcf_df['mitomap_gbcnt'] = [mitomap_poly_db.get(vt, ('0',))[0] for vt in var_tuples]
        vcf_df['mitomap_af'] = pd.to_numeric(vcf_df['mitomap_gbcnt'], errors='coerce') / 61134
        
        mitomap_disease_data = [mitomap_disease_db.get(vt, ('', '', '', '')) for vt in var_tuples]
        vcf_df['mitomap_status'] = [d[0] for d in mitomap_disease_data]
        vcf_df['mitomap_plasmy'] = [f"{d[1]}/{d[2]}" if d[1] or d[2] else '' for d in mitomap_disease_data]
        vcf_df['mitomap_disease'] = [d[3] for d in mitomap_disease_data]

        def get_haplo_status(row):
            variant_key = row['REF'] + row['POS'] + row['ALT']
            assoc_haplos = haplovar_db.get(variant_key)
            if assoc_haplos: return "haplo_var_match" if str(row["Haplogroup"]).lower() in assoc_haplos else "haplo_var_diff_haplo"
            return "not_haplo_var"
        vcf_df['Haplogroup_Var_Status'] = vcf_df.apply(get_haplo_status, axis=1)
        
        vcf_df['apogee_class'] = [str(apogee_db.get(vt, '')).strip("()").replace("'", "").replace(", ", "/") for vt in var_tuples]
        vcf_df['mitotip_class'] = [mitotip_db.get(vt, '') for vt in var_tuples]
        
        def get_hmtvar(vt):
            val = hmtvar_db.get(vt)
            if val and val[0]:
                try:
                    return json.loads(val[0]).get("pathogenicity", "")
                except json.JSONDecodeError:
                    return ""
            return ""
        vcf_df['hmtvar_class'] = [get_hmtvar(vt) for vt in var_tuples]
        
        all_variants_list.append(vcf_df)
    
    all_variants_df = pd.concat(all_variants_list, ignore_index=True)
    all_variants_df['AF'] = pd.to_numeric(all_variants_df['AF'], errors='coerce')
    
    print("[*] Counting variants per sample...")
    variant_counts = all_variants_df.groupby('SAMPLE_ID').apply(lambda df: pd.Series({'total_variants': len(df), 'heteroplasmic_variants': ((df['AF'] >= 0.05) & (df['AF'] <= 0.95)).sum(), 'homoplasmic_variants': (df['AF'] > 0.95).sum(), 'low_AF_variants': (df['AF'] < 0.05).sum(), 'haplogroup': df['Haplogroup'].iloc[0] if not df.empty else None})).reset_index()
    count_path = os.path.join(args.final_output_dir, f"sample_variant_count.txt")
    variant_counts.to_csv(count_path, sep='\t', index=False)
    print(f"[+] Variant count file saved to: {count_path}")

    print("[*] Performing final variant filtering...")
    
    all_variants_df['variant_key'] = all_variants_df['POS'].astype(str) + ':' + all_variants_df['REF'] + ':' + all_variants_df['ALT']

    # --- MODIFICATION: Calculate in_cohort_AC on ALL samples (minus exclusions) BEFORE filtering for POC ---
    print("[*] Calculating in-cohort allele count across all samples...")
    excluded_samples = ['370822', '378636', '378640']
    temp_df = all_variants_df[~all_variants_df['SAMPLE_ID'].str.startswith(tuple(excluded_samples))].copy()
    in_cohort_ac = temp_df['variant_key'].value_counts().reset_index()
    in_cohort_ac.columns = ['variant_key', 'in_cohort_AC']
    all_variants_df = all_variants_df.merge(in_cohort_ac, on='variant_key', how='left')
    all_variants_df['in_cohort_AC'] = all_variants_df['in_cohort_AC'].fillna(0).astype(int)
    
    print("[*] Adding maternal variant annotations...")
    proband_mother_map = {}
    for fam_id in master_pedigree_df['familyid'].unique():
        family_subset = master_pedigree_df[master_pedigree_df['familyid'] == fam_id]
        proband = family_subset[family_subset['category'] == 'POC']['sample_id']
        mother = family_subset[family_subset['category'] == 'Mother']['sample_id']
        if not proband.empty and not mother.empty:
            proband_full_id = all_variants_df[all_variants_df['SAMPLE_ID'].str.startswith(proband.iloc[0])]['SAMPLE_ID'].unique()
            mother_full_id = all_variants_df[all_variants_df['SAMPLE_ID'].str.startswith(mother.iloc[0])]['SAMPLE_ID'].unique()
            if len(proband_full_id) > 0 and len(mother_full_id) > 0:
                 proband_mother_map[proband_full_id[0]] = mother_full_id[0]

    all_variants_df['present_in_mother'] = False
    all_variants_df['mother_AF'] = pd.NA

    mother_variants_map = {mother_id: set(all_variants_df[all_variants_df['SAMPLE_ID'] == mother_id]['variant_key']) for mother_id in proband_mother_map.values()}
    mother_af_map = all_variants_df.set_index(['SAMPLE_ID', 'variant_key'])['AF'].to_dict()

    for proband_id, mother_id in proband_mother_map.items():
        proband_indices = all_variants_df[all_variants_df['SAMPLE_ID'] == proband_id].index
        mother_variants_set = mother_variants_map.get(mother_id, set())
        
        for idx in proband_indices:
            variant_key = all_variants_df.at[idx, 'variant_key']
            if variant_key in mother_variants_set:
                all_variants_df.at[idx, 'present_in_mother'] = True
                all_variants_df.at[idx, 'mother_AF'] = mother_af_map.get((mother_id, variant_key))

    poc_variants = all_variants_df[all_variants_df['Sample_Category'] == 'POC'].copy()
    
    # Calculate Freq (per-sample variant count) on the POC subset
    sample_counts = poc_variants['SAMPLE_ID'].value_counts().reset_index()
    sample_counts.columns = ['SAMPLE_ID', 'Freq']
    poc_variants = poc_variants.merge(sample_counts, on='SAMPLE_ID', how='left')

    final_column_order = [
        'SAMPLE_ID', 'variant_key', 'POS', 'REF', 'ALT', 'FILTER', 'GT', 'AD', 'AF', 'DP', 
        'Consequence', 'SYMBOL', 'BIOTYPE', 'HGVSc', 'HGVSp', 'Codons', 'VARIANT_CLASS', 
        'Haplogroup', 'Haplogroup_Var_Status', 'Sample_Category', 'gnomad_max_hl', 
        'gnomad_af_hom', 'gnomad_af_het', 'apogee_class', 'mitotip_class', 'hmtvar_class', 
        'helix_max_hl', 'helix_af_hom', 'helix_af_het', 'mitomap_gbcnt', 'mitomap_af', 
        'mitomap_status', 'mitomap_plasmy', 'mitomap_disease', 'clinvar_interp', 'in_cohort_AC', 
        'present_in_mother', 'mother_AF', 'Freq'
    ]
    for col in final_column_order:
        if col not in poc_variants.columns:
            poc_variants[col] = pd.NA
    poc_variants = poc_variants[final_column_order]

    prefilter_path = os.path.join(args.final_output_dir, f"POC_variant_list_prefiltering.txt")
    poc_variants.to_csv(prefilter_path, sep='\t', index=False, na_rep='')
    print(f"[+] Pre-filtering variant list saved to: {prefilter_path}")
    
    num_cols = ['gnomad_af_hom', 'helix_af_hom', 'mitomap_af', 'AF', 'Freq', 'in_cohort_AC']
    poc_variants[num_cols] = poc_variants[num_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    filtered_variants = poc_variants[
        (poc_variants['gnomad_af_hom'] < 0.01) &
        (poc_variants['helix_af_hom'] < 0.01) &
        (poc_variants['mitomap_af'] < 0.01) &
        (poc_variants['Consequence'] != 'synonymous_variant') &
        (poc_variants['Freq'] < 100) &
        (poc_variants['AF'] > 0.05) &
        (~poc_variants['clinvar_interp'].isin(['Benign', 'Likely benign'])) &
        (poc_variants['Haplogroup_Var_Status'] != 'haplo_var_match')
    ].copy()
    
    if not filtered_variants.empty:
        filtered_variants['mitomap_disease_sort'] = filtered_variants['mitomap_disease'].apply(lambda x: 0 if pd.isna(x) or x == '' else 1)
        filtered_variants = filtered_variants.sort_values(by=['mitomap_disease_sort', 'gnomad_af_hom', 'helix_af_hom', 'mitomap_af'], ascending=[False, True, True, True]).drop(columns=['mitomap_disease_sort'])

    if not filtered_variants.empty:
        filtered_variants = filtered_variants[final_column_order]

    final_path = os.path.join(args.final_output_dir, f"POC_variant_list.txt")
    filtered_variants.to_csv(final_path, sep='\t', index=False, na_rep='')
    print(f"[+] Final filtered variant list saved to: {final_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process and filter VEP annotated VCF files.")
    parser.add_argument("--vep-vcf-dir", required=True)
    parser.add_argument("--final-output-dir", required=True)
    parser.add_argument("--fullhaplogroups", required=True)
    parser.add_argument("--contamination", required=True)
    parser.add_argument("--master-pedigree-file", required=True)
    parser.add_argument("--gnomadcache", required=True)
    parser.add_argument("--clinvarcache", required=True)
    parser.add_argument("--mitomap-polycache", required=True)
    parser.add_argument("--mitomap-diseasecache", required=True)
    parser.add_argument("--helixcache", required=True)
    parser.add_argument("--haplogroup-varcache", required=True)
    parser.add_argument("--mitimpactcache", required=True)
    parser.add_argument("--mitotipcache", required=True)
    parser.add_argument("--hmtvarcache", required=True)
    
    args = parser.parse_args()
    process_and_filter_variants(args)

