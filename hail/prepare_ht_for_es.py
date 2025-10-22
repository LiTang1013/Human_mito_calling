#!/usr/bin/env python
import argparse
import hail as hl
import pprint

from utils import (
    get_expr_for_alt_allele,
    get_expr_for_contig,
    get_expr_for_consequence_lc_lof_flag,
    get_expr_for_variant_lc_lof_flag,
    get_expr_for_genes_with_lc_lof_flag,
    get_expr_for_consequence_loftee_flag_flag,
    get_expr_for_variant_loftee_flag_flag,
    get_expr_for_genes_with_loftee_flag_flag,
    get_expr_for_ref_allele,
    get_expr_for_variant_id,
    get_expr_for_vep_sorted_transcript_consequences_array,
    get_expr_for_xpos,
)


def reformat_vep_fields(ht):
    ht = ht.annotate(sortedTranscriptConsequences=get_expr_for_vep_sorted_transcript_consequences_array(vep_root=ht.vep))

    ht = ht.annotate(
        flags=hl.struct(
            lc_lof=get_expr_for_variant_lc_lof_flag(ht.sortedTranscriptConsequences),
            lof_flag=get_expr_for_variant_loftee_flag_flag(ht.sortedTranscriptConsequences),
            #lcr=ds.lcr,
            #segdup=ds.segdup,
        ),
        sortedTranscriptConsequences=hl.bind(
            lambda genes_with_lc_lof_flag, genes_with_loftee_flag_flag: ht.sortedTranscriptConsequences.map(
                lambda csq: csq.annotate(
                    flags=hl.struct(
                        lc_lof=get_expr_for_consequence_lc_lof_flag(csq),
                        lc_lof_in_gene=genes_with_lc_lof_flag.contains(csq.gene_id),
                        lof_flag=get_expr_for_consequence_loftee_flag_flag(csq),
                        lof_flag_in_gene=genes_with_loftee_flag_flag.contains(csq.gene_id),
                        nc_transcript=(csq.category == "lof") & (csq.lof == ""),
                    )
                )
            ),
            get_expr_for_genes_with_lc_lof_flag(ht.sortedTranscriptConsequences),
            get_expr_for_genes_with_loftee_flag_flag(ht.sortedTranscriptConsequences),
        ),
    )

    return ht


'''
def variant_id(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int = None):
    """
    Expression for computing <chrom>-<pos>-<ref>-<alt>. Assumes alleles were split.
    Args:
        max_length: (optional) length at which to truncate the <chrom>-<pos>-<ref>-<alt> string
    Return:
        string: "<chrom>-<pos>-<ref>-<alt>"
    """
    contig = normalized_contig(locus.contig)
    var_id = contig + "-" + hl.str(locus.position) + "-" + alleles[0] + "-" + alleles[1]

    if max_length is not None:
        return var_id[0:max_length]

    return var_id


def normalized_contig(contig: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.rbind(hl.str(contig).replace("^chr", ""), lambda c: hl.if_else(c == "MT", "M", c))

'''


def reformat_general_fields(ht):

    ht = ht.select_globals()


    # Derived top level fields
    ht = ht.annotate(
        alt=get_expr_for_alt_allele(ht),
        chrom=get_expr_for_contig(ht.locus),
        pos=ht.locus.position,
        ref=get_expr_for_ref_allele(ht),
        variant_id=get_expr_for_variant_id(ht),
        xpos=get_expr_for_xpos(ht.locus),
    )

    return ht


def reformat_freq_fields(ht):


    ht = ht.transmute(
    	ac=ht.AC,
        ac_het=ht.AC_het,
        ac_hom=ht.AC_hom,
        an=ht.AN,
        af=ht.AF,
        af_het=ht.AF_het,
        af_hom=ht.AF_hom,
        max_heteroplasmy=ht.max_hl
    )

    #ht = ht.transmute(**{f"{field}": expr_for_field_with_subpopulations(ht.info, field) for field in fields_per_subpopulation })

    return ht

'''
def prepare_ht_for_es(ds: hl.Table) -> hl.Table:

    ds = ds.select(
        # ID
        variant_id=variant_id(ds.locus, ds.alleles),
        reference_genome=ds.locus.dtype.reference_genome.name,
        chrom=normalized_contig(ds.locus.contig),
        pos=ds.locus.position,
        ref=ds.alleles[0],
        alt=ds.alleles[1],
        # Frequency
        an=ds.AN,
        ac_hom=ds.AC_hom,
        ac_het=ds.AC_het,
        max_heteroplasmy=ds.max_hl,
        vep=ds.vep
    )

    return ds
'''

def prepare_ht_for_es(ht):
    ht = reformat_general_fields(ht)
    ht = reformat_freq_fields(ht)
    ht = reformat_vep_fields(ht)
    
    ht = ht.expand_types().drop("locus", "alleles", "vep")
    #ht = ht.expand_types().drop("locus", "alleles")

    return ht


def main(args):

	ht = hl.read_table(args.ht)

	#pprint.pprint(ht.describe())
	#pprint.pprint(ht.show())


	ht = reformat_general_fields(ht)
	ht = reformat_freq_fields(ht)
	ht = reformat_vep_fields(ht)

	ht = ht.expand_types().drop("locus", "alleles", "vep")

	pprint.pprint(ht.describe())
	pprint.pprint(ht.show())

	ht.write(args.out,overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--ht', '--input', '-i', help='bgzipped VCF file (.vcf.bgz)', required=True)
    #parser.add_argument('--meta', '-m', help='Meta file containing sample population and sex', required=True)
    parser.add_argument('--out', '-o', help='Hail table output file name', required=True)

    args = parser.parse_args()
    main(args)

