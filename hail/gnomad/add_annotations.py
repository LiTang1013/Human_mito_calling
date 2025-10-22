#!/usr/bin/env python
import argparse
import hail as hl
import logging
import re
import sys

from collections import Counter
from os.path import dirname
from textwrap import dedent

from gnomad.utils.annotations import age_hists_expr
from gnomad.utils.reference_genome import add_reference_sequence
from gnomad.utils.vep import vep_struct_to_csq
from gnomad_qc.v3.resources.meta import meta
from gnomad.resources.grch38.gnomad import POPS
from gnomad.sample_qc.ancestry import POP_NAMES #check this
POP_NAMES["mid"] = "Middle Eastern"


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("add annotations")
logger.setLevel(logging.INFO)

def format_mt(mt_path: str, hl_threshold: float = 0.95, alt_threshold: float = 0.01) -> hl.MatrixTable:
    """Sets HL to zero if below heteroplasmy threshold, adds in genotype based on HL

    :param str mt_path: path to the MatrixTable
    :param float hl_threshold: heteroplasmy level threshold to determine homoplasmic variants
    :param float alt_threshold: heteroplasmy level threshold that must be reached to define the variant as an alternative allele
    :return: MatrixTable with GT
    :rtype: MatrixTable
    """

    # remove alt threshold
    logger.info('Reading in MT...')
    mt = hl.read_matrix_table(mt_path)

    mt = mt.rename({'VL': 'HL'})

    # TODO: rename artifact-prone-site filter in combine script
    # replace hyphens in filters with underscores
    mt = mt.annotate_rows(filters = mt.filters.map(lambda x: x.replace('-', '_')))

    # convert array of single string to actual array with by splitting string on semicolon
    mt = mt.annotate_entries(FT=hl.str(mt.FT)[2:-2].split(';'))

    # add back in GT based on hl_threshold
    mt = mt.annotate_entries(GT=(hl.case()
        .when((mt.HL < hl_threshold) & (mt.HL > 0.0), hl.parse_call("0/1"))
        .when(mt.HL >= hl_threshold, hl.parse_call("1/1"))
        .when(mt.HL == 0, hl.parse_call("0/0"))
        .default(hl.null(hl.tcall))))

    return mt


def add_variant_context(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """Adds variant context annotations to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :rtype: hl.MatrixTable
    """

    # read in variant context data
    vc_ht = hl.import_table("gs://gnomad-kristen/mitochondria/resources/variant_context/chrM_pos_ref_alt_context_categories.txt", impute=True)

    # split columns into separate annotations
    vc_ht = vc_ht.annotate(ref=vc_ht['POS.REF.ALT'].split('\.')[1],
                          alt=vc_ht['POS.REF.ALT'].split('\.')[2],
                          strand=vc_ht.Context_category.split('_')[-1],
                          variant=vc_ht.Context_category.split('_')[0])

    # rename and select certain columns
    vc_ht = vc_ht.rename({'MT_POS' : 'pos', 'Annotation' : 'region'})
    vc_ht = vc_ht.select('pos', 'ref', 'alt', 'strand', 'region', 'variant')

    # key by locus and allele
    vc_ht = vc_ht.key_by(locus=hl.locus("MT", vc_ht.pos, reference_genome='GRCh37'),
        alleles=[vc_ht.ref, vc_ht.alt])

    # annotate original mt with variant context information
    input_mt = input_mt.annotate_rows(**vc_ht[input_mt.locus, input_mt.alleles])
    input_mt = input_mt.annotate_rows(variant_type=hl.str(input_mt.variant) + "_" + hl.str(input_mt.strand))

    return input_mt


def add_gnomad_metadata(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """Adds gnomad metadata to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :rtype: hl.MatrixTable
    """

    # should version here be hard coded instead?
    genome_meta_ht = hl.read_table("gs://gnomad/metadata/genomes_v3.1/gnomad_v3.1_sample_qc_metadata.ht")

    genome_meta_struct = genome_meta_ht[input_mt.s]

    input_mt = input_mt.annotate_cols(release=genome_meta_struct.release,
                                      hard_filters=genome_meta_struct.sample_filters.hard_filters,
                                      research_project=genome_meta_struct.project_meta.research_project,
                                      project_id=genome_meta_struct.project_meta.project_id,
                                      product=genome_meta_struct.project_meta.product,
                                      sample_pi=genome_meta_struct.project_meta.sample_pi,
                                      sex_karyotype=genome_meta_struct.sex_imputation.sex_karyotype,
                                      age=hl.if_else(hl.is_defined(genome_meta_struct.project_meta.age), genome_meta_struct.project_meta.age, genome_meta_struct.project_meta.age_alt),
                                      broad_external=genome_meta_struct.project_meta.broad_external,
                                      pop=genome_meta_struct.population_inference.pop)

    return input_mt


def filter_by_copy_number(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """Calculates the mitochondrial copy number based on mean mitochondrial coverage and median nuclear coverage. Filters out samples with more extreme copy numbers.
    
    :param MatrixTable input_mt: MatrixTable
    :return: MatrixTable filtered to samples with a copy number of at least 50 and less than 500, number samples below 50 removed, number samples above 500 removed
    :rtype: hl.MatrixTable, int, int
    """

    # calculate mitochondrial copy number, if median autosomal coverage is not present default to 30x
    input_mt = input_mt.annotate_cols(mito_cn=2*input_mt.mt_mean_coverage/hl.if_else(hl.is_missing(input_mt.wgs_median_coverage), 30, input_mt.wgs_median_coverage))
    n_removed_below_cn = input_mt.aggregate_cols(hl.agg.count_where(input_mt.mito_cn < 50))
    n_removed_above_cn= input_mt.aggregate_cols(hl.agg.count_where(input_mt.mito_cn > 500))
    
    #logger.info(f'Removed {n_samples_removed} samples with MT copy number below 50 or greater than 500')

    # remove sample with a mitochondrial copy number below 50 or greater than 500
    input_mt = input_mt.filter_cols((input_mt.mito_cn >= 50) & (input_mt.mito_cn <= 500))
    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))
    
    return input_mt, n_removed_below_cn, n_removed_above_cn


def filter_by_contamination(input_mt: hl.MatrixTable, output_dir: str) -> hl.MatrixTable:
    """Calculates contamination and filters out samples with contamination above 2%
    
    :param MatrixTable input_mt: MatrixTable
    :param str output_dir: output directory to which results should be output

    :return: MatrixTable filtered to samples without contamination, number of contaminated samples removed
    :rtype: hl.MatrixTable, int
    """
    
    over_85_expr = (input_mt.HL >= 0.85) & (input_mt.FT == ["PASS"]) & input_mt.hap_defining_variant & ~hl.str(input_mt.filters).contains("artifact_prone_site")

    input_mt = input_mt.annotate_cols(over_85_mean = hl.agg.filter(over_85_expr, hl.agg.mean(input_mt.HL)),
                                     over_85_count = hl.agg.filter(over_85_expr, hl.agg.count_where(hl.is_defined(input_mt.HL))),
                                     bt_85_and_99_mean = hl.agg.filter(over_85_expr & (input_mt.HL <= 0.998), hl.agg.mean(input_mt.HL)),
                                     bt_85_and_99_count = hl.agg.filter(over_85_expr & (input_mt.HL <= 0.998), hl.agg.count_where(hl.is_defined(input_mt.HL))),
    )

    input_mt = input_mt.annotate_cols(contam_high_het = hl.if_else(input_mt.bt_85_and_99_count >= 3, 1 - input_mt.bt_85_and_99_mean, 1 - input_mt.over_85_mean))

    # if contam_high_het is nan, set to 0 (to avoid filtering out missing values which would be more common with haplogroups closer to reference haplogroup)
    input_mt = input_mt.annotate_cols(contam_high_het = hl.if_else(hl.is_nan(input_mt.contam_high_het), 0, input_mt.contam_high_het))

    # find samples on border of .02 that may flip between < 0.02 and > 0.02 from issues with floating point precision and maark these samples for removal
    epsilon = 0.000001
    border_samples = hl.literal(input_mt.aggregate_cols(hl.agg.filter((input_mt.contam_high_het > (.02 - epsilon)) & (input_mt.contam_high_het < (.02 + epsilon)), hl.agg.collect((input_mt.s)))))

    # save sample contamination information to separate file
    input_mt = input_mt.annotate_cols(keep=(input_mt.contamination < .02) & (input_mt.freemix_percentage < 2) & (input_mt.contam_high_het < .02) & ~border_samples.contains(input_mt.s))
    # save sample contamination information to separate file
    n_contaminated = input_mt.aggregate_cols(hl.agg.count_where(~input_mt.keep))
    logger.info(f'Removed {n_contaminated} samples with contamination above 2%')
    sample_data = input_mt.select_cols('contamination', 'freemix_percentage', 'contam_high_het', 'over_85_mean','over_85_count', 'bt_85_and_99_mean','bt_85_and_99_count', 'keep')
    data_export = sample_data.cols()
    data_export.export(f"{output_dir}/sample_contamination.tsv")

    # filter out filter with any contamination level above 2% heteroplasmy
    input_mt = input_mt.filter_cols(input_mt.keep)
    input_mt = input_mt.drop('keep')

    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))
    
    return input_mt, n_contaminated


def add_terra_metadata(input_mt: hl.MatrixTable, all_output: str) -> hl.MatrixTable:
    """Adds terra metadata to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :param str all_out: path to metadata file downloaded from terra

    :rtype: hl.MatrixTable
    """

    # add haplogroup and mutect/terra output annotations
    ht = hl.import_table(all_output, types={'contamination': hl.tfloat64, 'freemix_percentage': hl.tfloat64, 'mt_mean_coverage': hl.tfloat64, 'wgs_median_coverage': hl.tfloat64}, missing='').key_by('s')
    ht = ht.rename({'entity:participant_id': 'participant_id'})

    ht = ht.select("participant_id", "contamination", "freemix_percentage", "major_haplogroup", "wgs_median_coverage", "mt_mean_coverage")

    input_mt = input_mt.annotate_cols(**ht[input_mt.s])

    # annotate the high level haplogroup by taking the first letter with the exception of H and L haplogroups
    input_mt = input_mt.annotate_cols(hap=hl.if_else(input_mt.major_haplogroup.startswith("HV") | input_mt.major_haplogroup.startswith("L"),
        input_mt.major_haplogroup[0:2],
        input_mt.major_haplogroup[0]))

    return input_mt


def add_hap_defining(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """Adds bool on whether or not a variant is a haplogroup-defining variant to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :rtype: hl.MatrixTable
    """
    
    # TODO: move dataset location
    hap_defining_variants = hl.import_table("gs://gnomad-kristen/mitochondria/resources/phylotree/rCRS-centered_phylo_vars_final_update.txt")

    hap_defining = hl.literal(set(hap_defining_variants.variant.collect()))
    input_mt = input_mt.annotate_rows(variant_collapsed=input_mt.alleles[0] + hl.str(input_mt.locus.position) + input_mt.alleles[1])
    input_mt = input_mt.annotate_rows(hap_defining_variant=hap_defining.contains(input_mt.variant_collapsed))  # set hap_defining_variant to True or False

    return input_mt


def add_trna_predictions(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """Adds tRNA predictions on pathogenicity (PON-mt-tRNA and MitoTIP) to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :rtype: hl.MatrixTable
    """
    
    # add pon predictions
    pon_predictions = hl.import_table("gs://gnomad-kristen/mitochondria/resources/trna_predictions/pon_mt_trna_predictions_08_27_2020.txt")

    # if reference allele from fasta doesn't match Reference_nucleotide, pon trna is reporting the allele of opposite strand and need to get reverse complement for ref and alt
    add_reference_sequence(hl.get_reference('GRCh37'))
    pon_predictions = pon_predictions.annotate(ref=hl.get_sequence('MT', hl.int(pon_predictions.mtDNA_position), reference_genome='GRCh37'))
    pon_predictions = pon_predictions.annotate(alt=hl.if_else(pon_predictions.Reference_nucleotide == pon_predictions.ref, pon_predictions.New_nucleotide, hl.reverse_complement(pon_predictions.New_nucleotide)))
    pon_predictions = pon_predictions.key_by(variant_id=pon_predictions.ref + hl.str(pon_predictions.mtDNA_position) + pon_predictions.alt)
    input_mt = input_mt.annotate_rows(pon_mt_trna_prediction = pon_predictions[input_mt.variant_collapsed].Classification.lower().replace(' ', '_'),
        pon_ml_probability_of_pathogenicity = hl.float(pon_predictions[input_mt.variant_collapsed].ML_probability_of_pathogenicity))
    
    # add mitotip predictions
    mitotip_predictions = hl.import_table("gs://gnomad-kristen/mitochondria/resources/trna_predictions/mitotip_scores_08_27_2020.txt")
    mitotip_predictions = mitotip_predictions.key_by(variant_id=mitotip_predictions.rCRS + hl.str(mitotip_predictions.Position) + mitotip_predictions.Alt)
    input_mt = input_mt.annotate_rows(mitotip_score = hl.float(mitotip_predictions[input_mt.variant_collapsed].MitoTIP_Score))
    # set pathogenicity based on mitotip scores, classifications obtained from mitotip's website
    input_mt = input_mt.annotate_rows(mitotip_trna_prediction=(hl.case()
                                                      .when(input_mt.mitotip_score > 16.25, "likely_pathogenic")
                                                      .when((input_mt.mitotip_score <= 16.25) & (input_mt.mitotip_score > 12.66), "possibly_pathogenic")
                                                      .when((input_mt.mitotip_score <= 12.66) & (input_mt.mitotip_score >= 8.44), "possibly_benign")
                                                      .when((input_mt.mitotip_score < 8.44), "likely_benign")
                                                      .or_missing()))
    
    return input_mt


def get_indel_expr(input_mt: hl.MatrixTable) -> hl.expr.BooleanExpression:
    """Generates expression for filtering to indels that should be used to evaluate indel stacks
    
    :param MatrixTable input_mt: MatrixTable
    :return: indels to be used to evaluate indel stacks
    :rtype: hl.expr.BooleanExpression
    """
    
    indel_expr = hl.is_indel(input_mt.alleles[0], input_mt.alleles[1]) & (input_mt.HL <= 0.95) & (input_mt.HL >= 0.01) & (input_mt.FT == ["PASS"])

    return indel_expr 


def generate_expressions(input_mt: hl.MatrixTable, hl_threshold: float = 0.95) -> hl.MatrixTable:
    """Creates expressions to use for annotating the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :param float hl_threshold: heteroplasmy level threshold to determine homoplasmic variants
    :return: tuple of hail expressions
    :rtype: tuple
    """
 
    # calculate AC and AN
    AC = hl.agg.count_where((input_mt.HL > 0.0))
    AN = hl.agg.count_where(hl.is_defined(input_mt.HL))
    AF = AC/AN

    # calculate AC for het and hom variants, and histogram for HL
    AC_hom = hl.agg.count_where(input_mt.HL >= hl_threshold)
    AC_het = hl.agg.count_where((input_mt.HL < hl_threshold) & (input_mt.HL > 0.0))
    HL_hist = hl.agg.filter(input_mt.HL > 0, hl.agg.hist(input_mt.HL, 0, 1, 10))
    DP_hist_alt = hl.agg.filter(input_mt.GT.is_non_ref(), hl.agg.hist(input_mt.DP, 0, 2000, 10))
    DP_hist_all = hl.agg.hist(input_mt.DP, 0, 2000, 10)
    DP_mean = hl.agg.mean(input_mt.DP)
    MQ_mean = hl.agg.mean(input_mt.MQ)
    TLOD_mean = hl.agg.mean(input_mt.TLOD)

    # calculate AF
    AF_hom = AC_hom/AN
    AF_het = AC_het/AN

    # calculate max individual heteroplasmy
    max_HL = hl.agg.max(input_mt.HL)

    # haplogroup annotations
    pre_hap_AC = hl.agg.group_by(input_mt.hap, AC)
    pre_hap_AN = hl.agg.group_by(input_mt.hap, AN)
    pre_hap_AF = hl.agg.group_by(input_mt.hap, AF)
    pre_hap_AC_het = hl.agg.group_by(input_mt.hap, AC_het)
    pre_hap_AC_hom = hl.agg.group_by(input_mt.hap, AC_hom)
    pre_hap_AF_hom = hl.agg.group_by(input_mt.hap, AF_hom)
    pre_hap_AF_het = hl.agg.group_by(input_mt.hap, AF_het)
    pre_hap_HL_hist = hl.agg.group_by(input_mt.hap, HL_hist.bin_freq)
    pre_hap_FAF = hl.agg.group_by(input_mt.hap, hl.experimental.filtering_allele_frequency(hl.int32(AC), hl.int32(AN), 0.95))
    pre_hap_FAF_hom = hl.agg.group_by(input_mt.hap, hl.experimental.filtering_allele_frequency(hl.int32(AC_hom), hl.int32(AN), 0.95))

    # population annotations
    pre_pop_AC = hl.agg.group_by(input_mt.pop, AC)
    pre_pop_AN = hl.agg.group_by(input_mt.pop, AN)
    pre_pop_AF = hl.agg.group_by(input_mt.pop, AF)
    pre_pop_AC_het = hl.agg.group_by(input_mt.pop, AC_het)
    pre_pop_AC_hom = hl.agg.group_by(input_mt.pop, AC_hom)
    pre_pop_AF_hom = hl.agg.group_by(input_mt.pop, AF_hom)
    pre_pop_AF_het = hl.agg.group_by(input_mt.pop, AF_het)
    pre_pop_HL_hist = hl.agg.group_by(input_mt.pop, HL_hist.bin_freq)

    return hl.struct(AC=AC,
        AN=AN,
        AF=AF,
        AC_hom=AC_hom,
        AC_het=AC_het,
        hl_hist=HL_hist,
        dp_hist_all=DP_hist_all,
        dp_hist_alt=DP_hist_alt,
        dp_mean=DP_mean,
        mq_mean=MQ_mean,
        tlod_mean=TLOD_mean,
        AF_hom=AF_hom,
        AF_het=AF_het,
        max_hl=max_HL,
        pre_hap_AC=pre_hap_AC,
        pre_hap_AN=pre_hap_AN,
        pre_hap_AF=pre_hap_AF,
        pre_hap_AC_het=pre_hap_AC_het,
        pre_hap_AF_het=pre_hap_AF_het,
        pre_hap_AC_hom=pre_hap_AC_hom,
        pre_hap_AF_hom=pre_hap_AF_hom,
        pre_hap_hl_hist=pre_hap_HL_hist,
        pre_hap_faf=pre_hap_FAF,
        pre_hap_faf_hom=pre_hap_FAF_hom,
        pre_pop_AN=pre_pop_AN,
        pre_pop_AC_het=pre_pop_AC_het,
        pre_pop_AF_het=pre_pop_AF_het,
        pre_pop_AC_hom=pre_pop_AC_hom,
        pre_pop_AF_hom=pre_pop_AF_hom,
        pre_pop_hl_hist=pre_pop_HL_hist)


def standardize_haps(input_mt: hl.MatrixTable, annotation: str, haplogroup_order: list) -> list:
    """Converts the dictionary of haplogroup annotations into an array of values in a predefined haplogroup order

    :param MatrixTable input_mt: MatrixTable
    :param str annotation: annotation to convert and sort
    :param list haplogroup_order: order in which to sort the haplogroups
    :return: sorted list of haplogroup annotations (the values of the dictionary)
    :rtype: list
    """

    # converts haplogroup dictionary to sorted array
    value = [input_mt[annotation][x] for x in haplogroup_order]

    return value

def standardize_pops(input_mt: hl.MatrixTable, annotation: str, population_order: list) -> list:
    """Converts the dictionary of haplogroup annotations into an array of values in a predefined haplogroup order

    :param MatrixTable input_mt: MatrixTable
    :param str annotation: annotation to convert and sort
    :param list population_order: order in which to sort the populations
    :return: sorted list of population annotations (the values of the dictionary)
    :rtype: list
    """

    # converts haplogroup dictionary to sorted array
    value = [input_mt[annotation][x] for x in population_order]

    return value 


def add_variant_annotations(input_mt: hl.MatrixTable, hl_threshold: float = 0.95) -> hl.MatrixTable:
    """Adds variant annotations to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :param float hl_threshold: heteroplasmy level threshold to determine homoplasmic variants
    :return: annotated MatrixTable
    :rtype: hl.MatrixTable
    """

    # add variant annotations
    mt = input_mt.annotate_rows(**dict(generate_expressions(input_mt, hl_threshold)))

    # order the haplogroup-specific annotations
    list_hap_order = list(set(mt.hap.collect()))
    mt = mt.annotate_globals(hap_order=sorted(list_hap_order))

    # sanity check for haplogroups (make sure that they at least start with a letter) 
    for i in list_hap_order:
        if not re.match('^[A-Z]', i):
            sys.exit(f'Invalid haplogroup {i}, does not start with a letter')
        
    # generate histogram for site quality metrics across all variants
    # TODO: decide on bin edges
    dp_hist_all_variants = mt.aggregate_rows(hl.agg.hist(mt.dp_mean, 0, 4000, 40))
    mt = mt.annotate_globals(dp_hist_all_variants_bin_freq = dp_hist_all_variants.bin_freq,
                            dp_hist_all_variants_n_larger = dp_hist_all_variants.n_larger,
                            dp_hist_all_variants_bin_edges = dp_hist_all_variants.bin_edges)

    mq_hist_all_variants = mt.aggregate_rows(hl.agg.hist(mt.mq_mean, 0, 80, 40))  # is 80 the actual max value here?
    mt = mt.annotate_globals(mq_hist_all_variants_bin_freq = mq_hist_all_variants.bin_freq,
                            mq_hist_all_variants_n_larger = mq_hist_all_variants.n_larger,
                            mq_hist_all_variants_bin_edges = mq_hist_all_variants.bin_edges)
    
    tlod_hist_all_variants = mt.aggregate_rows(hl.agg.hist(mt.tlod_mean, 0, 40000, 40))  
    mt = mt.annotate_globals(tlod_hist_all_variants_bin_freq = tlod_hist_all_variants.bin_freq,
                            tlod_hist_all_variants_n_larger = tlod_hist_all_variants.n_larger,
                            tlod_hist_all_variants_bin_edges = tlod_hist_all_variants.bin_edges)

    # generate histogram for overall age distribution
    age_hist_all_samples = mt.aggregate_cols(hl.agg.hist(mt.age, 30, 80, 10)) 
    mt = mt.annotate_globals(age_hist_all_samples_bin_freq = age_hist_all_samples.bin_freq,
                            age_hist_all_samples_n_larger = age_hist_all_samples.n_larger,
                            age_hist_all_samples_n_smaller = age_hist_all_samples.n_smaller,
                            age_hist_all_samples_bin_edges = age_hist_all_samples.bin_edges)

    pre_hap_annotation_labels = ['pre_hap_AC',
                            'pre_hap_AN',
                            'pre_hap_AF',
                            'pre_hap_AC_het',
                            'pre_hap_AC_hom',
                            'pre_hap_AF_hom',
                            'pre_hap_AF_het',
                            'pre_hap_hl_hist',
                            'pre_hap_faf',
                            'pre_hap_faf_hom']

    for i in pre_hap_annotation_labels:
        final_annotation = re.sub("pre_", "", i)  # remove "pre" prefix for final annotations
        mt = mt.annotate_rows(**{final_annotation: standardize_haps(mt, i, sorted(list_hap_order))})

    # get a list of indexes where AC of the haplogroup is greater than 0, then get the list of haplogroups with that index
    mt = mt.annotate_rows(alt_haps=hl.zip_with_index(mt.hap_AC).filter(lambda x: x[1] > 0).map(lambda x: mt.hap_order[x[0]]))
    # count number of haplogroups containing an alt allele
    mt = mt.annotate_rows(n_alt_haps=hl.len(mt.alt_haps))

    # calculate hapmax
    mt = mt.annotate_rows(hapmax_AF_hom=mt.hap_order[(hl.argmax(mt.hap_AF_hom, unique=True))],
                          hapmax_AF_het=mt.hap_order[(hl.argmax(mt.hap_AF_het, unique=True))]
                         )

    # calculate faf hapmax
    mt = mt.annotate_rows(faf_hapmax=hl.max(mt.hap_faf),
                          faf_hapmax_hom=hl.max(mt.hap_faf_hom)
                         )
  
    # add populatation annotations
    found_pops = set(mt.pop.collect())
    final_pops = list(found_pops)
    final_pops = [x for x in POPS if x in final_pops]  # order according to POPS

    if len((found_pops) - set(POPS)) > 0:
        sys.exit(f'Invalid population found')
    mt = mt.annotate_globals(pop_order=final_pops)

    pre_pop_annotation_labels = ['pre_pop_AN',
                        'pre_pop_AC_het',
                        'pre_pop_AC_hom',
                        'pre_pop_AF_hom',
                        'pre_pop_AF_het',
                        'pre_pop_hl_hist']

    for i in pre_pop_annotation_labels:
        final_annotation = re.sub("pre_", "", i)  # remove "pre" prefix for final annotations
        mt = mt.annotate_rows(**{final_annotation: standardize_pops(mt, i, final_pops)})

    # add age histograms per variant type (heteroplasmic or homoplasmic)
    age_data = age_hists_expr(True, mt.GT, mt.age)
    mt = mt.annotate_rows(age_hist_hom=age_data.age_hist_hom,
                     age_hist_het=age_data.age_hist_het)
    
    # drop intermediate annotations
    annotations_to_drop = ["pre_hap_AC", 
                            "pre_hap_AN", 
                            "pre_hap_AF", 
                            "pre_hap_AC_het", 
                            "pre_hap_AC_hom", 
                            "pre_hap_AF_hom", 
                            "pre_hap_AF_het", 
                            "pre_hap_hl_hist", 
                            "pre_hap_faf", 
                            "pre_hap_faf_hom", 
                            "AC_mid_het", 
                            "AF_mid_het", 
                            "pre_pop_AN", 
                            "pre_pop_AC_het", 
                            "pre_pop_AC_hom", 
                            "pre_pop_AF_hom", 
                            "pre_pop_AF_het", 
                            "pre_pop_hl_hist"]

    mt = mt.drop(*annotations_to_drop)
    # last-minute drops (ever add back in?)
    mt = mt.drop("AC", "AF", "hap_AC", "hap_AF", "hap_faf", "faf_hapmax", "alt_haps", "n_alt_haps")

    mt = mt.annotate_rows(filters=hl.if_else(mt.filters=={"PASS"}, hl.empty_set(hl.tstr), mt.filters))

    return mt 


def generate_filter_histogram(input_mt: hl.MatrixTable, filter_name: str) -> hl.ArrayExpression:
    """Generates histogram for number of indiviudals with the specified sample-level filter at different heteroplasmy levels
    
    :param MatrixTable input_mt: MatrixTable
    :param str filter_name: name of sample-filter to count for histogram  
    :return: histogram containing the counts of individuals with a variant filtered by the specified filter name across binned heteroplasmy levels
    :rtype: hl.ArrayExpression
    """
    
    filter_histogram = hl.agg.filter(hl.str(input_mt.FT).contains(filter_name), hl.agg.hist(input_mt.HL, 0, 1, 10)).bin_freq
    
    return filter_histogram 


def add_filter_annotations(input_mt: hl.MatrixTable, alt_threshold: float = 0.01) -> hl.MatrixTable:
    """Generates histogram for number of individuals with the specified sample-level filter at different heteroplasmy levels
    
    :param MatrixTable input_mt: MatrixTable
    :param float alt_threshold: heteroplasmy level threshold that must be reached to define the variant as an alternative allele
    :return: MatrixTable with added annotations for sample and variant level filters and number of genotypes with heteroplasmy_below_10_percent
    :rtype: MatrixTable, int
    """

    # TODO: pull these from header instead?
    filters = ['base_qual',
               #'mt_many_low_hets',
               'position',
               #'possible_numt',
               'strand_bias',
               'weak_evidence',
               'contamination',
               'heteroplasmy_below_10_percent']

    # undo possible_numt and mt_many_low_hets filters
    filters_to_remove = hl.literal({"possible_numt", "mt_many_low_hets"})
    input_mt = input_mt.annotate_entries(FT=hl.array(hl.set(input_mt.FT).difference(filters_to_remove)))

    # if no filters exists after removing those specified above, set the FT field to PASS
    input_mt = input_mt.annotate_entries(FT=hl.if_else(hl.len(input_mt.FT) == 0, ["PASS"], input_mt.FT))

    # add common_low_heteroplasmy flag to variants where the overall frequency is > 0.001 for samples with a heteroplasmy > 0 and < 0.50 and either "low_allele_frac" or "PASS" for the genotype filter
    input_mt = input_mt.checkpoint("gs://gnomad-kristen/mitochondria/gnomadv3_1/gnomad_10_22_2020/test_laf.mt", overwrite=True)  # full matrix table for internal use
    input_mt = input_mt.annotate_rows(AC_mid_het=hl.agg.count_where((input_mt.HL < 0.50) & (input_mt.HL > 0.0) & ((input_mt.FT == ["PASS"]) | (input_mt.FT == ["low_allele_frac"]))))
    input_mt = input_mt.annotate_rows(AF_mid_het=input_mt.AC_mid_het/hl.agg.count_where(hl.is_defined(input_mt.HL) & ((input_mt.FT == ["PASS"]) | (input_mt.FT == ["low_allele_frac"]))))
    input_mt = input_mt.annotate_rows(common_low_heteroplasmy=input_mt.AF_mid_het > 0.001)
    common_low_het_count = input_mt.aggregate_rows(hl.agg.count_where(input_mt.common_low_heteroplasmy))
    logger.info(f'The "common_low_heteroplasmy" flag was applied to {common_low_het_count} variants')

    # set HL to 0 if < alt_threshold and remove variants that no longer have at least one alt call
    input_mt = input_mt.annotate_entries(HL=hl.if_else((input_mt.HL > 0) & (input_mt.HL < alt_threshold), 0, input_mt.HL))
    # filter field for all variants with a heteroplasmy of 0 should be set to PASS
    # this step is needed to prevent homref calls that are filtered
    input_mt = input_mt.annotate_entries(FT=hl.if_else(input_mt.HL < alt_threshold, ["PASS"], input_mt.FT))
    input_mt = input_mt.annotate_entries(GT=hl.if_else(input_mt.HL < alt_threshold, hl.parse_call("0/0"), input_mt.GT))

    # check that variants no longer contain the "low_allele_frac" filter (alt_threshold should be set to appropriate level to remove these variants)
    laf_rows = input_mt.filter_rows(hl.agg.any(hl.str(input_mt.FT).contains("low_allele_frac")))
    n_laf_rows = laf_rows.count_rows()
    if n_laf_rows > 0:
        sys.exit("low_allele_frac filter should no longer be present after applying alt_threshold (alt_threshold should equal the vaf_cutoff set in the mutect wdl)")
    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))

    # add variant-level indel_stack at any indel allele where all samples with a variant call had at least 2 different indels called at that position
    # if any sample had a solo indel at that position, do not filter 
    indel_expr = get_indel_expr(input_mt)
    input_mt = input_mt.annotate_cols(indel_pos_counter=hl.agg.filter(indel_expr, hl.agg.counter(input_mt.locus.position)))
    input_mt = input_mt.annotate_cols(indel_pos_counter_keys=input_mt.indel_pos_counter.keys())
    indel_expr = get_indel_expr(input_mt)
    input_mt = input_mt.annotate_entries(indel_occurences=(hl.case()
                                                                .when((indel_expr & (input_mt.indel_pos_counter.get(input_mt.locus.position) >= 2)), "stack")
                                                                .when((indel_expr & (input_mt.indel_pos_counter.get(input_mt.locus.position) == 1)), "solo")
                                                                .or_missing()))

    #if stack true and solo false, stack only and the indel should be filtered out
    input_mt = input_mt.annotate_rows(filters=hl.if_else(hl.agg.any(input_mt.indel_occurences=="stack") & ~hl.agg.any(input_mt.indel_occurences=="solo"), input_mt.filters.add("indel_stack"), input_mt.filters))
    indel_stack_count = input_mt.aggregate_rows(hl.agg.count_where(hl.str(input_mt.filters).contains("indel_stack")))
    logger.info(f'The "indel_stack" filter was applied to {indel_stack_count} variants')

    # add sample level filter to remove variants below 10% heteroplasmy
    input_mt = input_mt.annotate_entries(FT=hl.if_else((input_mt.HL < 0.10) & (input_mt.GT.is_het()), input_mt.FT.append("heteroplasmy_below_10_percent"), input_mt.FT))
    n_het_below_10_perc = input_mt.aggregate_entries(hl.agg.count_where(hl.str(input_mt.FT).contains("heteroplasmy_below_10_percent")))
    logger.info(f'The "heteroplasmy_below_10_percent" filter was applied to {n_het_below_10_perc} genotypes')

    # add npg (no pass genotypes) filter to sites that don't have at least one pass alt call 
    input_mt = input_mt.annotate_rows(filters=hl.if_else(~(hl.agg.any((input_mt.HL > 0.0) & (input_mt.FT == ["PASS"]))), input_mt.filters.add("npg"), input_mt.filters))
    npg_count = input_mt.aggregate_rows(hl.agg.count_where(hl.str(input_mt.filters).contains("npg")))
    logger.info(f'The "npg" filter was applied to {npg_count} variants')

    for i in filters:
        annotation_name = i + "_hist"
        input_mt = input_mt.annotate_rows(**{annotation_name: generate_filter_histogram(input_mt, i)})
    input_mt = input_mt.annotate_rows(excluded_AC=hl.agg.count_where(input_mt.FT != ["PASS"]))

    # remove "PASS" from filters column if it's not the only filter
    input_mt = input_mt.annotate_rows(filters=hl.if_else(input_mt.filters != {"PASS"}, input_mt.filters.remove("PASS"), input_mt.filters))
    
    return input_mt, n_het_below_10_perc


def filter_genotypes(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """Sets all genotype field values to missing if the variant is not "PASS" for that sample
    
    :param MatrixTable input_mt: MatrixTable
    :return: MatrixTable with filtered genotype fields set to missing
    :rtype: MatrixTable
    """

    pass_expr = (input_mt.FT == ["PASS"])
    
    input_mt = input_mt.annotate_entries(GT = hl.or_missing(pass_expr, input_mt.GT),
                             DP = hl.or_missing(pass_expr, input_mt.DP),
                             HL = hl.or_missing(pass_expr, input_mt.HL),
                             FT = hl.or_missing(pass_expr, input_mt.FT),
                             MQ = hl.or_missing(pass_expr, input_mt.MQ),
                             TLOD = hl.or_missing(pass_expr, input_mt.TLOD))
    
    return input_mt


def add_sample_annotations(input_mt: hl.MatrixTable, hl_threshold: float = 0.95) -> hl.MatrixTable:
    """Adds sample annotations to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :param float hl_threshold: heteroplasmy level threshold to determine homoplasmic variants 
    :return: MatrixTable with sample annotations
    :rtype: hl.MatrixTable
    """

    # count number of variants
    num_rows = input_mt.count_rows()

    # add sample qc annotations
    no_artifact_expr = ~hl.str(input_mt.filters).contains("artifact_prone_site")                             
    input_mt = input_mt.annotate_cols(callrate=hl.agg.filter(no_artifact_expr, (hl.agg.count_where(hl.is_defined(input_mt.HL)))/num_rows),
                         n_singletons_het=hl.agg.filter(no_artifact_expr, hl.agg.count_where((input_mt.AC_het == 1) & ((input_mt.HL < hl_threshold) & (input_mt.HL > 0.0)))),
                         n_singletons_hom=hl.agg.filter(no_artifact_expr, hl.agg.count_where((input_mt.AC_hom == 1) & (input_mt.HL >= hl_threshold))),
                         n_snp_het=hl.agg.filter(no_artifact_expr, hl.agg.count_where((input_mt.HL < hl_threshold) & (input_mt.HL > 0.0) &  hl.is_snp(input_mt.alleles[0], input_mt.alleles[1]))),
                         n_snp_hom=hl.agg.filter(no_artifact_expr, hl.agg.count_where((input_mt.HL >= hl_threshold) & hl.is_snp(input_mt.alleles[0], input_mt.alleles[1]))),
                         n_indel_het=hl.agg.filter(no_artifact_expr, hl.agg.count_where((input_mt.HL < hl_threshold) & (input_mt.HL > 0.0) & (~hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])))),
                         n_indel_hom=hl.agg.filter(no_artifact_expr, hl.agg.count_where((input_mt.HL >= hl_threshold) & (~hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])))))
    
    return input_mt


def add_vep(input_mt: hl.MatrixTable, run_vep: bool, vep_output: str) -> hl.MatrixTable:
    """Adds vep annotations to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :param bool run_vep: whether or not to run vep
    :param str vep_output: path to the MatrixTable output vep results (either the existing results or where to ouput new vep results)
    :return: MatrixTable with vep annotations
    :rtype: hl.MatrixTable
    """

    if run_vep:
        vep_mt = hl.vep(input_mt)
        vep_mt = vep_mt.checkpoint(vep_output, overwrite=True)
    else:
        vep_mt = hl.read_matrix_table(vep_output)

    input_mt = input_mt.annotate_rows(vep=vep_mt.index_rows(input_mt.locus, input_mt.alleles).vep)
    # TODO: get vep version directly from config file
    input_mt = input_mt.annotate_globals(vep_version="v101")


    # if only filter is END_TRUNC, change lof for LC to HC and remove the END_TRUNC filter
    # remove SINGLE_EXON flags because all exons are single exon in the mitochondria
    input_mt = input_mt.annotate_rows(
        vep = input_mt.vep.annotate(
            transcript_consequences = 
            input_mt.vep.transcript_consequences.map(lambda x: x.annotate(
                lof=hl.if_else(
                    x.lof_filter == "END_TRUNC", "HC", x.lof),
                lof_filter = hl.if_else(
                    x.lof_filter == "END_TRUNC", hl.null(hl.tstr), x.lof_filter),
                lof_flags=hl.if_else(
                    x.lof_flags == "SINGLE_EXON", hl.null(hl.tstr), x.lof_flags)))))

    end_trunc_count = input_mt.filter_rows(hl.str(input_mt.vep.transcript_consequences[0].lof_filter).contains("END_TRUNC")).count_rows()
    if end_trunc_count > 0:
        sys.exit(f"END_TRUNC filter should no longer be present but was found for {end_trunc_count} variants")

    single_exon_count = input_mt.filter_rows(hl.str(input_mt.vep.transcript_consequences[0].lof_flags).contains("SINGLE_EXON")).count_rows()
    if single_exon_count > 0:
        sys.exit(f"SINGLE_EXON flag should no longer be present but was found for {single_exon_count} variants")

    return input_mt

def add_rsids(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """Adds rsid annotations to the MatrixTable

    :param MatrixTable input_mt: MatrixTable
    :return: MatrixTable with rsid annotations
    :rtype: hl.MatrixTable
    """
    
    DBSNP_B154_CHR_CONTIG_RECODING = {'NC_012920.1': 'chrM'}

    import_args={
        "path": "gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.bgz",
        "header_file": "gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.header",
        "force_bgz": True,
        "contig_recoding": DBSNP_B154_CHR_CONTIG_RECODING,
        "skip_invalid_loci": True,
        "min_partitions": 400,
        "reference_genome": "GRCh38"}

    dbsnp = hl.import_vcf(**import_args)
    input_mt = input_mt.annotate_rows(rsid=dbsnp.index_rows(input_mt.locus, input_mt.alleles).rsid)
    input_mt = input_mt.annotate_globals(dbsnp_version="b154")

    
    return input_mt


def export_simplified_variants(input_ht: hl.Table,  output_dir: str) -> None:
    """Exports only several high-level variant annotations
    
    :param hl.Table input_ht: Hail Table of variants
    :param str output_dir: output directory to which results should be output
    :return: None
    """
    
    reduced_ht = input_ht.key_by(chromosome = input_ht.locus.contig,
                                 position = input_ht.locus.position,
                                 ref = input_ht.alleles[0],
                                 alt = input_ht.alleles[1]).select('filters', 'AC_hom', 'AC_het', 'AF_hom', 'AF_het', 'AN', 'max_hl').rename({'max_hl' : 'max_observed_heteroplasmy'})
    reduced_ht = reduced_ht.annotate(filters=hl.if_else(hl.len(reduced_ht.filters)==0, "PASS", hl.str(',').join(hl.array(reduced_ht.filters))))

    reduced_ht.export(f'{output_dir}/reduced_annotations.txt')

def report_stats(input_mt: hl.MatrixTable, output_dir: str, pass_only: bool, n_samples_below_cn: int, n_samples_above_cn: int, n_samples_contam: int, n_het_below_10_perc: int, min_threshold: float = 0.10, hl_threshold: float = 0.95) -> None:
    """Generates output report with basic stats
    
    :param MatrixTable input_mt: MatrixTable
    :param str output_dir: output directory to which results should be output
    :param bool pass_only: whether or not directory should be filtered to pass_only variants
    :param int n_samples_below_cn: number of samples removed because mitochondrial number is less than 50
    :param int n_samples_above_cn: number of samples removed because mitochondrial number is above 500
    :param int n_samples_contam: number of samples removed because of contamination
    :param int n_het_below_10_perc: number of genotypes filtered because the heteroplasmy levels was below 10%
    :param float min_threshold: minimum threshold to pass genotypes
    :param float hl_threshold: heteroplasmy level threshold to determine homoplasmic variants 
    :return: None
    """
    
    if pass_only:
        suffix = "_pass"
        input_mt = input_mt.filter_rows(hl.len(input_mt.filters) == 0)
    else:
        suffix = ""
    out_stats = hl.hadoop_open(f'{output_dir}/stats{suffix}.txt', 'w')

    if pass_only:
        out_stats.write("Below metrics are for PASS-only variants\n\n")
    
    # report numbers of filtered samples/genotypes
    out_stats.write(f'Number of samples removed because contamination above 2%: {n_samples_contam}\n')
    out_stats.write(f'Number of samples removed because mitochondrial copy number below 50: {n_samples_below_cn}\n')
    out_stats.write(f'Number of samples removed because mitochondrial copy number above 500: {n_samples_above_cn}\n')
    out_stats.write(f'Number of genotypes filtered because "heteroplasmy_below_10_percent": {n_het_below_10_perc}\n\n')

    # count variant, samples, bases
    unique_variants, samples = input_mt.count()
    out_stats.write(f'Number of samples: {samples}\n')
    out_stats.write(f'Number of unique variants: {unique_variants}\n')
    bases_w_variant = len(set(input_mt.locus.position.collect()))
    out_stats.write(f'Number of bases with variation: {bases_w_variant}\n\n')

    # count number of filters
    for filter_name, filter_count in Counter([i for sublist in input_mt.filters.collect() for i in sublist]).items():
        out_stats.write(f'Number of variants with "{filter_name}" filter: {filter_count} variants\n')

    
    # calculate row stats
    row_stats = input_mt.aggregate_rows(hl.struct(common_low_het_count=hl.agg.count_where(input_mt.common_low_heteroplasmy),
        het_only_sites=hl.agg.count_where((input_mt.AC_het > 0) & (input_mt.AC_hom == 0)),
        snps=hl.agg.count_where(hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])), 
        indels=hl.agg.count_where(hl.is_indel(input_mt.alleles[0], input_mt.alleles[1]))))

    # calculate entry stats
    entry_stats=input_mt.aggregate_entries(hl.struct(total_variants=hl.agg.count_where(input_mt.HL > 0),
                                        total_hom_variants=hl.agg.count_where(input_mt.HL >= hl_threshold),
                                        total_het_variants=hl.agg.count_where((input_mt.HL < hl_threshold) & (input_mt.HL >= min_threshold)),
                                        min_hl=hl.agg.filter(input_mt.HL > 0, hl.agg.min(input_mt.HL)),
                                        max_hl=hl.agg.filter(input_mt.HL > 0, hl.agg.max(input_mt.HL))))

    # count number of flags
    out_stats.write(f'Number of variants with "common_low_heteroplasmy" flag: {row_stats["common_low_het_count"]} variants\n\n')

    # count variants
    out_stats.write(f'Total number of variants: {entry_stats["total_variants"]}\n')

    # count of homoplasmic and heteroplasmic variants
    out_stats.write(f'Number of heteroplasmic-only sites: {row_stats["het_only_sites"]}\n')

    percent_hom = round(entry_stats["total_hom_variants"]/entry_stats["total_variants"], 2)
    percent_het = round(entry_stats["total_het_variants"]/entry_stats["total_variants"], 2)
    out_stats.write(f'Total number of homoplasmic variants: {entry_stats["total_hom_variants"]}\n')
    out_stats.write(f'Percent homoplasmic variants: {percent_hom}\n')
    out_stats.write(f'Total number of heteroplasmic variants: {entry_stats["total_het_variants"]}\n')
    out_stats.write(f'Percent heteroplasmic variants: {percent_het}\n\n')

    out_stats.write(f'Minimum heteroplasmy detected: {entry_stats["min_hl"]}\n')
    out_stats.write(f'Maximum heteroplasmy detected: {entry_stats["max_hl"]}\n\n')

    # count number of snps and indels
    out_stats.write(f'Number of SNPs: {row_stats["snps"]}\n')
    out_stats.write(f'Number of indels: {row_stats["indels"]}\n')
    
    out_stats.close()

def change_to_grch38_chrm(input_mt: hl.MatrixTable) -> None:
    """Changes build to GRCh38 and filters reference genome to chrM

    :param MatrixTable input_mt: MatrixTable
    :return: MatrixTable with GRCh38 reference genomed subsetted to just chrM
    :rtype: hl.MatrixTable
    """

    ref = hl.get_reference("GRCh38")
    my_ref = hl.ReferenceGenome("GRCh38_chrM", contigs=['chrM'],  lengths={'chrM': ref.lengths['chrM']})
    assert('chrM' in ref.contigs)
    input_mt = input_mt.key_rows_by(locus=hl.locus("chrM", input_mt.locus.position, reference_genome='GRCh38_chrM'), alleles=input_mt.alleles)
    
    return input_mt

def format_vcf(input_mt: hl.MatrixTable, output_dir: str, hl_threshold: float = 0.95, alt_threshold: float = 0.01) -> dict:
    """Generate dictionary for vcf header annotations

    :param MatrixTable input_mt: MatrixTable
    :param float hl_threshold: heteroplasmy level threshold to determine homoplasmic variants
    :param float alt_threshold: heteroplasmy level threshold that must be reached to define the variant as an alternative allele
    :param str output_dir: output directory to which appended header info should be written
    :return: MatrixTable with vcf annotations in the info field and dictionary of filter, info, and format fields to be output in the vcf header; path of vcf headers to append
    :rtype: hl.MatrixTable, dict, str
    """
    input_mt = change_to_grch38_chrm(input_mt)

    haplogroup_order = hl.eval(input_mt.hap_order)
    population_order = hl.eval(input_mt.pop_order)

    age_hist_hom_bin_edges = input_mt.age_hist_hom.bin_edges.take(1)[0]
    age_hist_het_bin_edges = input_mt.age_hist_het.bin_edges.take(1)[0]
    hl_hist_bin_edges = input_mt.hl_hist.bin_edges.take(1)[0]
    dp_hist_all_bin_edges = input_mt.dp_hist_all.bin_edges.take(1)[0]
    dp_hist_alt_bin_edges = input_mt.dp_hist_alt.bin_edges.take(1)[0]

    input_mt = input_mt.annotate_rows(hl_hist=input_mt.hl_hist.bin_freq,
                    age_hist_hom_bin_freq=input_mt.age_hist_hom.bin_freq,
                    age_hist_hom_n_smaller=input_mt.age_hist_hom.n_smaller,
                    age_hist_hom_n_larger=input_mt.age_hist_hom.n_larger,
                    age_hist_het_bin_freq=input_mt.age_hist_het.bin_freq,
                    age_hist_het_n_smaller=input_mt.age_hist_het.n_smaller,
                    age_hist_het_n_larger=input_mt.age_hist_het.n_larger,
                    dp_hist_all_n_larger=input_mt.dp_hist_all.n_larger,
                    dp_hist_alt_n_larger=input_mt.dp_hist_alt.n_larger,
                    dp_hist_all_bin_freq=input_mt.dp_hist_all.bin_freq,
                    dp_hist_alt_bin_freq=input_mt.dp_hist_alt.bin_freq)
    
    input_mt = input_mt.annotate_rows(vep=vep_struct_to_csq(input_mt.vep))

    # get length of annotations to use in Number fields in the vcf where necessary
    len_hap_hl_hist = len(input_mt.hap_hl_hist.take(1)[0])
    len_pop_hl_hist = len(input_mt.pop_hl_hist.take(1)[0])

    # output appended header info to file
    vcf_header_file = output_dir + "/extra_fields_for_header.tsv"
    appended_vcf_header = dedent(f"""
    ##VEP version={hl.eval(input_mt.vep_version)}
    ##dbSNP version={hl.eval(input_mt.dbsnp_version)}
    ##age distributions=bin_edges:{hl.eval(input_mt.age_hist_all_samples_bin_edges)}, bin_freq:{hl.eval(input_mt.age_hist_all_samples_bin_freq)}, n_smaller:{hl.eval(input_mt.age_hist_all_samples_n_smaller)}, n_larger:{hl.eval(input_mt.age_hist_all_samples_n_larger)}

    """)
    with hl.hadoop_open(vcf_header_file, 'w') as out:
        out.write(appended_vcf_header)
 
    # drop intermediate annotations
    input_mt = input_mt.drop('qual', 'pos', 'ref', 'alt', 'strand', 'region', 'variant', 'variant_type', 'age_hist_het', 'age_hist_hom', 'dp_hist_all', 'dp_hist_alt')

    ht = input_mt.rows()
    input_mt = input_mt.annotate_rows(info=input_mt.info.annotate(**ht[input_mt.row_key]).drop('info', 'AF', 'rsid')).select_rows("rsid", "filters", "info")  # create info annotation (and drop nested info annotation)

    # convert "," to "|" for array annotations
    for key,value in input_mt.row_value.info.items():
        if str(value).startswith("<Array"):
            if str(value.dtype).startswith("array<array"):
                # if value is an array of arrays, only replace the commas within each individual array
                input_mt = input_mt.annotate_rows(info=input_mt.info.annotate(**{key: hl.map(lambda x: hl.delimit(x, delimiter="|"), input_mt.info[key])}))
            else:
                input_mt = input_mt.annotate_rows(info=input_mt.info.annotate(**{key: hl.delimit(input_mt.info[key], delimiter="|")}))

    meta_dict = {'filter': {
                    'artifact_prone_site': {'Description': 'Variant overlaps site that is commonly reported in literature to be artifact prone'},
                    'npg': {'Description': 'No-pass-genotypes site (no individuals were PASS for the variant)'},
                    'indel_stack': {'Description': 'Allele where all samples with the variant call had at least 2 different indels called at the position'}
    },
                 'info': {
                     'variant_collapsed': {'Description': 'Variant in format of RefPosAlt', 'Number': '1', 'Type': 'String'}, 
                     'hap_defining_variant': {'Description': 'Present if variant is present as a haplogroup defining variant in PhyloTree build 17', 'Number': '0', 'Type': 'Flag'},
                     'common_low_heteroplasmy': {'Description': f'Present if variant is found at an overall frequency of .001 across all samples with a heteroplasmy level > 0 and < 0.50 (includes variants <{alt_threshold} heteroplasmy which are subsequently filtered)', 'Number': '0', 'Type': 'Flag'},
                     'AN': {'Description': 'Overall allele number (number of samples with non-missing genotype)', 'Number': '1', 'Type': 'Integer'},
                     'AC_hom': {'Description': f'Allele count restricted to variants with a heteroplasmy level >= {hl_threshold}', 'Number': '1', 'Type': 'Integer'},  # should put in threshold variable    
                     'AC_het': {'Description': f'Allele count restricted to variants with a heteroplasmy level >= 0.10 and < {hl_threshold}', 'Number': '1', 'Type': 'Integer'},
                     'hl_hist': {'Description': f'Histogram of heteroplasmy levels; bin edges are: {hl_hist_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'hap_hl_hist': {'Description': f'Histogram of heteroplasmy levels for each haplogroup; bin edges are: {hl_hist_bin_edges}, haplogroup order: {haplogroup_order}', 'Number': f'{len_hap_hl_hist}', 'Type': 'String'},
                     'AF_hom': {'Description': f'Allele frequency restricted to variants with a heteroplasmy level >= {hl_threshold}', 'Number': '1', 'Type': 'Float'},
                     'AF_het': {'Description': f'Allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < {hl_threshold}', 'Number': '1', 'Type': 'Float'},
                     'max_hl': {'Description': 'Maximum heteroplasmy level observed among all samples for that variant', 'Number': '1', 'Type': 'Float'},
                     'hap_AN': {'Description': f'List of overall allele number for each haplogroup, haplogroup order: {haplogroup_order}', 'Number': '1', 'Type': 'String'},
                     'hap_AC_het': {'Description': f'List of AC_het for each haplogroup, haplogroup order: {haplogroup_order}', 'Number': '1', 'Type': 'String'},
                     'hap_AC_hom': {'Description': f'List of AC_hom for each haplogroup, haplogroup order: {haplogroup_order}', 'Number': '1', 'Type': 'String'},
                     'hap_AF_hom': {'Description': f'List of AF_hom for each haplogroup, haplogroup order: {haplogroup_order}', 'Number': '1', 'Type': 'String'},
                     'hap_AF_het': {'Description': f'List of AF_het for each haplogroup, haplogroup order: {haplogroup_order}', 'Number': '1', 'Type': 'String'},
                     'hap_faf_hom': {'Description': f'List of filtering allele frequency for each haplogroup restricted to homoplasmic variants, haplogroup order: {haplogroup_order}', 'Number': '1', 'Type': 'String'},
                     'hapmax_AF_hom': {'Description': 'Haplogroup with maximum AF_hom', 'Number': '1', 'Type': 'String'},
                     'hapmax_AF_het': {'Description': 'Haplogroup with maximum AF_het', 'Number': '1', 'Type': 'String'},
                     'faf_hapmax_hom': {'Description': 'Maximum filtering allele frequency across haplogroups restricted to homoplasmic variants', 'Number': '1', 'Type': 'Float'},
                     'bin_edges_hl_hist': {'Description': 'Bin edges for histogram of heteroplasmy levels', 'Number': '1', 'Type': 'String'},
                     'pop_AN': {'Description': f'List of overall allele number for each population, population order: {population_order}', 'Number': '1', 'Type': 'String'},
                     'pop_AC_het': {'Description': f'List of AC_het for each population, population order: {population_order}', 'Number': '1', 'Type': 'String'},
                     'pop_AC_hom': {'Description': f'List of AC_hom for each population, population order: {population_order}', 'Number': '1', 'Type': 'String'},
                     'pop_AF_hom': {'Description': f'List of AF_hom for each population, population order: {population_order}', 'Number': '1', 'Type': 'String'},
                     'pop_AF_het': {'Description': f'List of AF_het for each population, population order: {population_order}', 'Number': '1', 'Type': 'String'},
                     'pop_hl_hist': {'Description': f'Histogram of heteroplasmy levels for each population; bin edges are: {hl_hist_bin_edges}, population order: {population_order}', 'Number': f'{len_pop_hl_hist}', 'Type': 'String'},
                     'age_hist_hom_bin_freq': {'Description': f'Histogram of ages of individuals with a homoplasmic variant; bin edges are: {age_hist_hom_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'age_hist_hom_n_smaller': {'Description': 'Count of age values falling below lowest histogram bin edge for individuals with a homoplasmic variant', 'Number': '1', 'Type': 'String'},
                     'age_hist_hom_n_larger': {'Description': 'Count of age values falling above highest histogram bin edge for individuals with a homoplasmic variant', 'Number': '1', 'Type': 'String'},
                     'age_hist_het_bin_freq': {'Description': f'Histogram of ages of individuals with a heteroplasmic variant; bin edges are: {age_hist_het_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'age_hist_het_n_smaller': {'Description': 'Count of age values falling below lowest histogram bin edge for individuals with a heteroplasmic variant', 'Number': '1', 'Type': 'String'},
                     'age_hist_het_n_larger': {'Description': 'Count of age values falling above highest histogram bin edge for individuals with a heteroplasmic variant', 'Number': '1', 'Type': 'String'},
                     'dp_hist_all_n_larger': {'Description': 'Count of dp values falling above highest histogram bin edge for all individuals', 'Number': '1', 'Type': 'String'},
                     'dp_hist_alt_n_larger': {'Description': 'Count of dp values falling above highest histogram bin edge for individuals with the alternative allele', 'Number': '1', 'Type': 'String'},
                     'dp_hist_all_bin_freq': {'Description': f'Histogram of dp values for all individuals; bin edges are: {dp_hist_all_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'dp_hist_alt_bin_freq': {'Description': f'Histogram of dp values for individuals with the alternative allele; bin edges are: {dp_hist_alt_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'dp_mean': {'Description': 'Mean depth across all individuals for the site', 'Number': '1', 'Type': 'Float'},
                     'mq_mean': {'Description': 'Mean MMQ (median mapping quality) across individuals with a variant for the site', 'Number': '1', 'Type': 'Float'},
                     'tlod_mean': {'Description': 'Mean TLOD (Log 10 likelihood ratio score of variant existing versus not existing) across individuals with a variant for the site', 'Number': '1', 'Type': 'Float'},
                     'pon_mt_trna_prediction': {'Description': 'tRNA pathogenicity classification from PON-mt-tRNA', 'Number': '1', 'Type': 'String'},
                     'pon_ml_probability_of_pathogenicity': {'Description': 'tRNA ML_probability_of_pathogenicity from PON-mt-tRNA', 'Number': '1', 'Type': 'Float'},
                     'mitotip_score': {'Description': 'MitoTip raw score', 'Number': '1', 'Type': 'Float'},
                     'mitotip_trna_prediction': {'Description': 'MitoTip score interpretation', 'Number': '1', 'Type': 'String'},
                     'vep': {'Description': 'Consequence annotations from Ensembl VEP; note that the SINGLE_EXON flag and END_TRUNC filters have been removed from the LOFTEE annotations to avoid misinterpretation in context of the mitochondrial genome. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info', 'Number': '.', 'Type': 'String'},
                     'filters': {'Description': 'Site-level filters', 'Number': '.', 'Type': 'String'},
                     'base_qual_hist': {'Description': f'Histogram of number of individuals failing the base_qual filter (alternate allele median base quality) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'heteroplasmy_below_10_percent_hist': {'Description': f'Histogram of number of individuals with a heteroplasmy level below 10%, bin edges are: {hl_hist_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'position_hist': {'Description': f'Histogram of number of individuals failing the position filter (median distance of alternate variants from end of reads) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'strand_bias_hist': {'Description': f'Histogram of number of individuals failing the strand_bias filter (evidence for alternate allele comes from one read direction only) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'weak_evidence_hist': {'Description': f'Histogram of number of individuals failing the weak_evidence filter (mutation does not meet likelihood threshold) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'contamination_hist': {'Description': f'Histogram of number of individuals failing the contamination filter across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}', 'Number': '1', 'Type': 'String'},
                     'excluded_AC': {'Description': 'Excluded allele count (number of individuals in which the variant was filtered out)', 'Number': '1', 'Type': 'String'}
                },
             'format': {
                 'GT': {'Description': f'Genotype, 1/1 if heteroplasmy level >= {hl_threshold}, and 0/1 if heteroplasmy level < {hl_threshold}', 'Number': '1', 'Type': 'String'},
                 'DP': {'Description': 'Depth of coverage', 'Number': '1', 'Type': 'Integer'},
                 'FT': {'Description': 'Sample-level filters', 'Number': '.', 'Type': 'String'},
                 'HL': {'Description': 'Heteroplasmy level', 'Number': '1', 'Type': 'Float'},
                 'MQ': {'Description': 'Mapping quality', 'Number': '1', 'Type': 'Float'},
                 'TLOD': {'Description': 'Log 10 likelihood ratio score of variant existing versus not existing', 'Number': '1', 'Type': 'Float'}
                        }
            }
    
    return input_mt, meta_dict, vcf_header_file


def main(args):
    mt_path = args.mt_path
    all_output = args.all_output
    run_vep = args.run_vep
    vep_results = args.vep_results
    gnomad_subset = args.subset_to_gnomad_release
    hl_threshold = args.hl_threshold
    alt_threshold = args.alt_threshold

    logger.info(f'Cutoff for homoplasmic variants is set to {hl_threshold}...')
    
    # define mt path, output directory, subset name
    output_dir = dirname(mt_path)
    subset_name = ""

    logger.info('Reformatting MT...')
    mt = format_mt(mt_path, hl_threshold, alt_threshold)

    logger.info('Adding annotations from terra...')
    mt = add_terra_metadata(mt, all_output)

    logger.info('Annotating haplogroup-defining variants...')
    mt = add_hap_defining(mt)

    logger.info('Annotating tRNA predictions...')
    mt = add_trna_predictions(mt)

    logger.info('Adding gnomAD metadata sample annotations...')
    mt = add_gnomad_metadata(mt)

    logger.info('Adding variant context annotations...')
    mt = add_variant_context(mt)

    # if specified, subet to only the gnomAD samples in the current release
    if gnomad_subset:
        logger.warn('Subsetting results to gnomAD release samples...')
        subset_name = "_gnomad"

        # subset to release samples and filter out rows that no longer have at least one alt call
        mt = mt.filter_cols(mt.release)  # filter to cols where release is true
        mt = mt.filter_rows(hl.agg.any(mt.HL > 0))

    logger.info('Filtering out low copy number samples...')
    mt, n_removed_below_cn, n_removed_above_cn = filter_by_copy_number(mt)

    logger.info('Filtering out contaminated samples...')
    mt, n_contaminated = filter_by_contamination(mt, output_dir)

    logger.info('Switch build and checkpoint...')    
    # switch build 37 to build 38
    mt = mt.key_rows_by(locus=hl.locus("chrM", mt.locus.position, reference_genome='GRCh38'), alleles=mt.alleles)
    mt = mt.checkpoint(f"{output_dir}/prior_to_vep.mt", overwrite=args.overwrite)

    logger.info('Adding vep annotations...')
    mt = add_vep(mt, run_vep, vep_results)

    logger.info('Adding dbsnp annotations...')
    mt = add_rsids(mt)

    # set up output paths for callset
    annotated_mt_path = f'{output_dir}/annotated_combined{subset_name}.mt'
    sites_ht_path = f'{output_dir}/combined_sites_only{subset_name}.ht'   
    sites_txt_path = f'{output_dir}/combined_sites_only{subset_name}.txt'
    sites_vcf_path = f'{output_dir}/combined_sites_only{subset_name}.vcf.bgz'
    sample_txt_path = f'{output_dir}/sample_annotations{subset_name}.txt'
    samples_vcf_path = f'{output_dir}/sample_annotations{subset_name}.vcf.bgz'

    # output paths for callset with haplogroup-defining variants dropped
    annotated_mt_drop_hap_path = f'{output_dir}/annotated_combined_drop_hap{subset_name}.mt'
    sites_ht_drop_hap_path = f'{output_dir}/combined_sites_only_drop_hap{subset_name}.ht'
    sites_txt_drop_hap_path = f'{output_dir}/combined_sites_only_drop_hap{subset_name}.txt'
    sample_txt_drop_hap_path = f'{output_dir}/sample_annotations_drop_hap{subset_name}.txt'
    samples_vcf_drop_hap_path = f'{output_dir}/sample_annotations_drop_hap{subset_name}.vcf.bgz'

    logger.info('Results will be output to the following files:')
    print(annotated_mt_path)
    print(sites_ht_path)
    print(sites_txt_path)
    print(annotated_mt_drop_hap_path)
    print(sites_ht_drop_hap_path)
    print(sites_txt_drop_hap_path)
    print(sample_txt_path)
    print(sample_txt_drop_hap_path)
    print(samples_vcf_path)
    print(sites_vcf_path)
    print(samples_vcf_drop_hap_path)

    logger.info('Annotating MT...')
    mt, n_het_below_10_perc = add_filter_annotations(mt, alt_threshold)
    mt = filter_genotypes(mt)
    mt = add_variant_annotations(mt, hl_threshold)
    mt = mt.checkpoint(annotated_mt_path, overwrite=args.overwrite)  # full matrix table for internal use

    logger.info('Generating summary statistics reports...')
    report_stats(mt, output_dir, False, n_removed_below_cn, n_removed_above_cn, n_contaminated, n_het_below_10_perc)
    report_stats(mt, output_dir, True, n_removed_below_cn, n_removed_above_cn, n_contaminated, n_het_below_10_perc)

    variant_ht = mt.rows()
    variant_ht = variant_ht.drop('ref', 'alt', 'pos', 'strand', 'region', 'variant', 'variant_type', 'info')
    variant_ht.export(sites_txt_path)  # sites-only txt file for external use
    variant_ht.write(sites_ht_path, overwrite=args.overwrite)  # sites-only ht for external use
    mt = add_sample_annotations(mt, hl_threshold)

    sample_ht = mt.cols()
    sample_ht.group_by(sample_ht.hap).aggregate(n=hl.agg.count()).export(f'{output_dir}/haplogroup_counts.txt')  # counts of top level haplogroups
    sample_ht.export(sample_txt_path)  # sample annotations txt file for internal use
    
    rows_ht = mt.rows()
    export_simplified_variants(rows_ht, output_dir)
    vcf_mt, vcf_meta, vcf_header_file = format_vcf(mt, output_dir, hl_threshold)
    hl.export_vcf(vcf_mt, samples_vcf_path, metadata=vcf_meta, append_to_header=vcf_header_file, tabix=True)  # full vcf for internal use
    vcf_variant_ht = vcf_mt.rows()
    rows_mt = hl.MatrixTable.from_rows_table(vcf_variant_ht).key_cols_by(s='foo')
    hl.export_vcf(rows_mt, sites_vcf_path, metadata=vcf_meta, append_to_header=vcf_header_file, tabix=True)  # sites-only vcf for external use

    if args.drop_hap_defining:
        logger.info('Annotating MT after dropping haplogroup-defining variants...')
        mt_drop_hap = mt.filter_rows(~mt.hap_defining_variant)  # remove haplogroup-defining variants
     
        mt_drop_hap = add_variant_annotations(mt_drop_hap, hl_threshold)
        mt_drop_hap, n_het_below_10_perc = add_filter_annotations(mt_drop_hap)
        mt_drop_hap = mt_drop_hap.checkpoint(annotated_mt_drop_hap_path, overwrite=args.overwrite)  # full matrix table for internal use (dropped haplogroup-defining variants)

        variant_ht = mt_drop_hap.rows()
        variant_ht.export(sites_txt_drop_hap_path)  # sites-only txt for internal use (dropped haplogroup-defining variants)
        variant_ht.write(sites_ht_drop_hap_path, overwrite=args.overwrite)  # sites-only txt for internal use (dropped haplogroup-defining variants)
        mt_drop_hap = add_sample_annotations(mt_drop_hap, hl_threshold) 

        sample_ht = mt_drop_hap.cols()
        sample_ht.export(sample_txt_drop_hap_path)  # sample annotations txt file for internal use (dropped haplogroup-defining variants)
        
        vcf_mt, vcf_meta, vcf_header_file = format_vcf(mt_drop_hap, output_dir, hl_threshold)
        hl.export_vcf(vcf_mt, samples_vcf_drop_hap_path, metadata=vcf_meta, append_to_header=vcf_header_file)  # full vcf for internal use (dropped haplogroup-defining variants)

    logger.info('All annotation steps are completed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script adds variant annotations to the mitochondria VCF/MT')
    parser.add_argument('-a', '--all_output', help='Output file that results from terra data download')
    parser.add_argument('-p', '--mt_path', help='Path to combined mt')
    parser.add_argument('-r', '--vep_results', help='MatrixTable path to output vep results (either the existing results or where to ouput new vep results)')
    parser.add_argument('-t', '--hl_threshold', help='Cutoff to define a variant as homoplasmic', nargs='?', const=1, type=float, default=0.95)
    parser.add_argument('-m', '--alt_threshold', help='Minimum cutoff to define a variant as an alternative allele', nargs='?', const=1, type=float, default=0.01)
    parser.add_argument('--subset_to_gnomad_release', help='Set to True to only include released gnomAD samples', action='store_true')
    parser.add_argument('--run_vep', help='Set to True to run/rerun vep', action='store_true')
    parser.add_argument('--overwrite', help='Overwrites existing files', action='store_true')
    parser.add_argument('--drop_hap_defining', help='Set to True to output file without haplogroup-defining variants', action='store_true')


    args = parser.parse_args()
    
    main(args)


