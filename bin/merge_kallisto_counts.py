#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path
import gffpandas.gffpandas as gffpd
from collections import Counter


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metadata_f",
        help = "TSV file containing metadata"
    )
    parser.add_argument(
        "--gff_f",
        help = "GFF annotation file"
    )
    parser.add_argument("--out_dir", help="Directory for results")
    args = parser.parse_args()
    merge_counts(**vars(args))

def merge_counts(metadata_f, gff_f, out_dir):
    # merge counts
    metadata = pd.read_csv(metadata_f, sep = "\t")
    quant_dfs = []
    for index, row in metadata.iterrows():
        sample_name = row['sample']
        quant_file = os.path.join('kallisto_'+sample_name, 'abundance.tsv')
        quant_dat = pd.read_csv(quant_file, sep = "\t")
        quant_dat = quant_dat[["target_id", "est_counts"]]
        quant_dat = quant_dat.rename(columns={'target_id': sample_name})
        quant_dfs.append(quant_dat)
    quant_dfs = [df.set_index('Gene') for df in quant_dfs]
    quant_merged = pd.concat(quant_dfs, axis=1)
    quant_merged.rename_axis('feature_id').reset_index()
    outf1 = os.path.join(out_dir, 'kallisto_merged_counts.tsv')
    quant_merged.to_csv(outf1, index=False, sep='\t')
    # make gene annotation df
    annot_dat = (gffpd.read_gff3(gff_f)).df
    cds_annot = annot_dat[annot_dat['type']=='CDS']
    annot_split = cds_annot['attributes'].str.split(";")
    locus_tags = (annot_split.str[1]).str.replace(r'Parent=gene-','')
    gene_names = ((annot_split.str[6]).str.split("=")).str[1]
    cds_annot['locus_tag'] = locus_tags
    cds_annot['gene_name'] = gene_names
    cds_annot['gene_length'] = (cds_annot['end'] - cds_annot['start']) + 1
    cds_annot['biotype'] = 'protein_coding'
    cds_annot = cds_annot[['locus_tag', 'biotype', 'gene_name', 'gene_length']]
    outf2 = os.path.join(out_dir, 'ref_gene_df.tsv')
    cds_annot.to_csv(outf2, index=False, sep='\t')


if __name__ == "__main__":
    parse()
