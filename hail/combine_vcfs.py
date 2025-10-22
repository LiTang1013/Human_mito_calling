import argparse
import logging
import math
import os
import re
import sys

import hail as hl
import pandas as pd

#from gnomad.resources.resource_utils import DataException
from hail.utils.java import info
from subprocess import Popen, PIPE, check_output
from typing import Dict

META_DICT = {
    "filter": {
        "artifact_prone_site": {
            "Description": "Variant overlaps an artifact_prone_site"
        }
    },
    "format": {
        "DP": {"Description": "Depth of coverage", "Number": "1", "Type": "Integer"},
        "FT": {"Description": "Sample-level filters", "Number": ".", "Type": "String"},
        "HL": {"Description": "Heteroplasmy level", "Number": "1", "Type": "Float"},
        "MQ": {"Description": "Mapping quality", "Number": "1", "Type": "Float"},
        "TLOD": {
            "Description": "Log 10 likelihood ratio score of variant existing versus not existing",
            "Number": "1",
            "Type": "Float",
        },
    },
}


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("combine_mitochondria_vcfs_into_mt")
logger.setLevel(logging.INFO)


def check_vcf_existence(
    participant_data: str, vcf_col: str, sample_map: str, output_bucket: str
) -> Dict[str, str]:
    """For each participant specified in sample_map, checks that the vcf file exists, and if so, add the sample and vcf path to a dictionary

    :param str participant_data: participant data (downloaded data tab from terra)
    :param str vcf_col: name of column that contains vcf output
    :param str sample_map: path to file of samples to subset (tab-delimited participant_id and sample)
    :param str output_bucket: path to bucket to which results should be written

    :return: dictionary of samples for which the vcf existence was confirmed (sample as key, path to vcf as value)
    :rtype: Dict[str, str]
    """

    # create file that will contain the samples with confirmed vcfs and their paths
    out_vcf = hl.hadoop_open(f"{output_bucket}/vcfs_to_combine.list", "w")

    # create participants_of_interest dictionary which will contain samples to which the results shoudl be subset
    participants_of_interest = {}
    confirmed_vcfs = {}
    with hl.hadoop_open(sample_map, "r") as f:
        next(f)
        for line in f:
            line = line.rstrip()
            items = line.split("\t")
            participant, sample = items[0:2]
            participants_of_interest[participant] = 0

    # load in data from terra
    participant_info = hl.import_table(participant_data)
    df = participant_info.to_pandas()

    # check if the sample is in participants_of_interest, check that the vcf exists, and if yes to both, add to confirmed_vcfs dictionary
    for _, row in df.iterrows():
        participant_id = row["entity:participant_id"]
        sample = row["s"]
        vcf = row[vcf_col]

        if participant_id in participants_of_interest and vcf != "":
            if hl.hadoop_is_file(vcf):
                out_vcf.write(f"{sample}\t{vcf}\n")
                confirmed_vcfs[sample] = vcf

    out_vcf.close()

    return confirmed_vcfs


def multi_way_union_mts(mts: list, tmp_dir: str, chunk_size: int) -> hl.MatrixTable:
    """Joins MatrixTables in the provided list

    :param list mts: list of MatrixTables to join together
    :param str tmp_dir: path to temporary directory for intermediate results
    :param int chunk_size: number of MatrixTables to join per chunk

    :return: joined MatrixTable
    :rtype: MatrixTable
    """

    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []
        for i in range(n_jobs):
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )

            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            next_stage.append(
                merged.checkpoint(
                    os.path.join(tmp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"done stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


def join_mitochondria_vcfs_into_mt(
    input_tsv: str, output_bucket: str, chunk_size: int = 100
) -> hl.MatrixTable:
    """Reformats and joins individual mitochondrial vcfs into one MatrixTable

    :param dict confirmed_vcfs: dictionary of samples for which the vcf existence was confirmed (sample as key, path to vcf as value)
    :param str output_bucket: path to bucket to which results should be written
    :param int chunk_size: number of MatrixTables to join per chunk

    :return: joined MatrixTable of samples given in confirmed_vcfs dictionary
    :rtype: hl.MatrixTable
    """

    mt_list = []

    with open(input_tsv, "r") as f:

        for line in f:
            line = line.rstrip()
            items = line.split("\t")
            sample, vcf_path = items[0:2]

            mt = hl.import_vcf(vcf_path, reference_genome="GRCh38")
            # because the vcfs are split, there is only one AF value, although misinterpreted as an array because Number=A in vcf header
            # second value of MMQ is the value for the alternate allele
            mt = mt.select_entries("DP", HL=mt.AF[0])
            mt = mt.annotate_entries(
                MQ=hl.float(mt.info["MMQ"][1]), TLOD=mt.info["TLOD"][0], FT=hl.if_else(hl.len(mt.filters) == 0, {"PASS"}, mt.filters)
            )
            # use GRCh37 as reference as this is more compatibile with mitochondria resources that may be added as annotations in downstream scripts
            mt = mt.key_rows_by(
                locus=hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
                alleles=mt.alleles,
            )
            mt = mt.key_cols_by(s=sample)
            mt = mt.select_rows()
            mt_list.append(mt)


    temp_out_dir = output_bucket + "/temp"
    combined_mt = multi_way_union_mts(mt_list, temp_out_dir, chunk_size)

    return combined_mt

def remove_FT_values(mt: hl.MatrixTable, filters_to_remove: list = ['possible_numt', 'mt_many_low_hets', 'FAIL', 'blacklisted_site']) -> hl.MatrixTable:
    """Removes the FT filters specified in filters_to_remove
    
    By default, this function removes the 'possible_numt', 'mt_many_low_hets', and 'FAIL' filters (because these filters were found to have low performance), 
    and the 'blacklisted_site' filter because this filter did not always behave as expected in early GATK versions (can be replaced with apply_mito_artifact_filter function)

    :param hl.MatrixTable mt:  MatrixTable
    :param list filters_to_remove: list of FT filters that should be removed from the entries
    
    :return: MatrixTable with certain FT filters removed
    :rtype: MatrixTable
    """
    
    filters_to_remove = hl.set(filters_to_remove)
    mt = mt.annotate_entries(FT=hl.array((mt.FT).difference(filters_to_remove)))

    # if no filters exists after removing those specified above, set the FT field to PASS
    mt = mt.annotate_entries(FT=hl.if_else(hl.len(mt.FT) == 0, ["PASS"], mt.FT))
    
    return(mt)

def determine_hom_refs(
    mt: hl.MatrixTable, coverage_mt_path: str, minimum_homref_coverage: int = 100
) -> hl.MatrixTable:
    """Uses coverage to distinguish between homref and missing sites, outputs the resulting mt

    :param hl.MatrixTable mt: MatrixTable to use an input
    :param str coverage_mt_path: MatrixTable of sample level coverage at each position
    :param int minimum_homref_coverage: minimum depth of coverage required to call a genotype homoplasmic reference rather than missing

    :return: MatrixTable with missing genotypes converted to homref depending on coverage
    :rtype: hl.MatrixTable
    """

    # convert coverage to build GRCh37
    # note: the mitochondrial reference genome is the same for GRCh38 and GRCh37
    coverages = hl.read_matrix_table(coverage_mt_path)
    coverages = coverages.key_rows_by(
        locus=hl.locus("MT", coverages.locus.position, reference_genome="GRCh37")
    )

    mt = mt.annotate_entries(
        DP=hl.if_else(hl.is_missing(mt.HL), coverages[mt.locus, mt.s].coverage, mt.DP)
    )

    hom_ref_expr = hl.is_missing(mt.HL) & (mt.DP > minimum_homref_coverage)

    mt = mt.annotate_entries(
        HL=hl.if_else(hom_ref_expr, 0.0, mt.HL),
        FT=hl.if_else(hom_ref_expr, ["PASS"], mt.FT),
        DP=hl.if_else(
            hl.is_missing(mt.HL) & (mt.DP <= minimum_homref_coverage),
            hl.null(hl.tint32),
            mt.DP,
        ),
    )

    return mt


def apply_mito_artifact_filter(
    mt: hl.MatrixTable,
    artifact_prone_sites_path: str,
) -> hl.MatrixTable:
    """Add back in artifact_prone_site filter

    :param hl.MatrixTable mt: MatrixTable to use an input
    :param str artifact_prone_sites_path: path to BED file of artifact_prone_sites to flag in the filters column

    :return: MatrixTable with artifact_prone_sites filter
    :rtype: hl.MatrixTable

    """

    # apply "artifact_prone_site" filter to any SNP or deletion that spans a known problematic site
    mt = mt.annotate_rows(
        position_range=hl.range(
            mt.locus.position, mt.locus.position + hl.len(mt.alleles[0])
        )
    )

    artifact_sites = []

    with open(artifact_prone_sites_path, "r") as f:
        for line in f:
            pos = line.split()[2]
            artifact_sites.append(int(pos))
    sites = hl.literal(set(artifact_sites))

    mt = mt.annotate_rows(
        filters=hl.if_else(
            hl.len(hl.set(mt.position_range).intersection(sites)) > 0,
            {"artifact_prone_site"},
            {"PASS"},
        )
    )

    mt = mt.drop("position_range")

    return mt


def main(args):

    '''
    participant_data = args.participant_data
    sample_map = args.sample_map
    coverage_mt_path = args.coverage_mt_path
    vcf_col = args.vcf_col
    artifact_prone_sites_path = args.artifact_prone_sites_path
    output_bucket = args.output_bucket
    file_suffix = args.file_suffix
    minimum_homref_coverage = args.minimum_homref_coverage
    chunk_size = args.chunk_size
    overwrite = args.overwrite
    '''

    input_tsv = args.input_tsv

    artifact_prone_sites_path = args.artifact_prone_sites_path
    chunk_size = args.chunk_size
    overwrite = args.overwrite
    output_bucket = args.output_bucket
    file_suffix = args.file_suffix

    minimum_homref_coverage = args.minimum_homref_coverage
    coverage_mt_path = args.coverage_mt_path

    '''
    logger.info("Confirming existence of individual sample vcfs...")
    confirmed_vcfs = check_vcf_existence(
        participant_data, vcf_col, sample_map, output_bucket
    )
    '''

    logger.info("Combining VCFs...")
    combined_mt = join_mitochondria_vcfs_into_mt(
        input_tsv, output_bucket, chunk_size
    )
    output_path_mt = f"{output_bucket}/raw_combined_mt.mt"
    combined_mt = combined_mt.checkpoint(output_path_mt, overwrite=overwrite)

    
    logger.info("Removing certain FT filters...")
    combined_mt = remove_FT_values(combined_mt)

    logger.info("Determining homoplasmic reference sites...")
    combined_mt = determine_hom_refs(
        combined_mt, coverage_mt_path, minimum_homref_coverage
    )
    
    logger.info("Applying artifact_prone_site fiter...")
    combined_mt = apply_mito_artifact_filter(combined_mt, artifact_prone_sites_path)

    logger.info("Writing combined MT and VCF...")
    # set the file names for output files
    out_vcf = f"{output_bucket}/combined_{file_suffix}.vcf.bgz"
    out_mt = f"{output_bucket}/combined_{file_suffix}.mt"

    combined_mt.write(out_mt, overwrite=True)
    hl.export_vcf(combined_mt, out_vcf, metadata=META_DICT)
    
    logger.info("DONE")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="This script combines individual mitochondria vcf files into one MatrixTable, determines homoplasmic reference sites, and applies an artifact_prone_site filter"
    )
    '''
    p.add_argument(
        "-p",
        "--participant_data",
        help="Participant data (downloaded data tab from terra)",
    )
    p.add_argument(
        "-s",
        "--sample_map",
        help="Path to file of samples to subset (tab-delimited participant_id and sample)",
    )
    '''

    p.add_argument(
        "-c", "--coverage_mt_path", help="Path to MatrixTable of sample-level coverage"
    )
    p.add_argument("-v", "--vcf_col", help="Name of column that contains vcf output")
    p.add_argument(
        "-a",
        "--artifact_prone_sites_path",
        help="List of artifact_prone sites to flag in the FILTER column (in BED format)",
    )


    p.add_argument(
        "-i",
        "--input_tsv",
        help="Input tsv file",
    )

    p.add_argument(
        "-o",
        "--output_bucket",
        help="Path to bucket to which results should be written",
    )

    p.add_argument(
        "-f",
        "--file_suffix",
        help="File suffix to append to names of output files",
    )
    p.add_argument(
        "--minimum_homref_coverage",
        help="Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing",
        nargs="?",
        const=1,
        type=int,
        default=100,
    )
    p.add_argument(
        "--chunk_size",
        help="Chunk size to use for combining vcfs",
        nargs="?",
        const=1,
        type=int,
        default=100,
    )
    p.add_argument("--overwrite", help="Overwrites existing files", action="store_true")

    args = p.parse_args()

    main(args)
