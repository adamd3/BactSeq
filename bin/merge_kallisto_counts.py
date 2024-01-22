#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path

# import gffpandas.gffpandas as gffpd
from collections import Counter


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_f", help="TSV file containing metadata")
    # parser.add_argument("--gff_f", help="GFF annotation file")
    parser.add_argument("--out_dir", help="Directory for results")
    args = parser.parse_args()
    merge_counts(**vars(args))


def merge_counts(
    metadata_f,
    # gff_f,
    out_dir,
):
    # merge counts
    metadata = pd.read_csv(metadata_f, sep="\t")
    quant_dfs = []
    for index, row in metadata.iterrows():
        sample_name = row["sample"]
        quant_file = os.path.join("kallisto_" + sample_name, "abundance.tsv")
        quant_dat = pd.read_csv(quant_file, sep="\t")
        quant_dat = quant_dat[["target_id", "est_counts", "length"]]
        quant_dat = quant_dat.rename(columns={"est_counts": sample_name})
        quant_dfs.append(quant_dat)
    quant_dfs = [df.set_index("target_id") for df in quant_dfs]
    quant_merged = pd.concat(quant_dfs, axis=1)
    quant_merged = quant_merged.loc[:, ~quant_merged.columns.duplicated()]
    # get gene lengths from kallisto results (ensures no divergence from GFF)
    quant_merged = quant_merged.rename_axis("feature_id").reset_index()
    gene_lengths = quant_merged[["feature_id", "length"]]
    gene_lengths = gene_lengths.rename(columns={"feature_id": "locus_tag"})
    quant_merged = quant_merged[["feature_id"] + metadata["sample"].tolist()]
    # export merged counts
    outf1 = os.path.join(out_dir, "gene_counts1.tsv")
    quant_merged.to_csv(outf1, index=False, sep="\t")
    # # gene annotation df for extracting protein-coding genes
    # annot_dat = (gffpd.read_gff3(gff_f)).df
    # cds_annot = (annot_dat[annot_dat['type']=='CDS']).copy(deep=True)
    # # rRNA_annot = annot_dat[annot_dat['type']=='rRNA']
    # annot_split = cds_annot['attributes'].str.split(";")
    # locus_tags = (annot_split.str[1]).str.replace(r'Parent=gene-','')
    # gene_names = ((annot_split.str[6]).str.split("=")).str[1]
    # cds_annot['locus_tag'] = locus_tags
    # cds_annot['gene_name'] = gene_names
    # cds_annot['gene_length'] = (cds_annot['end'] - cds_annot['start']) + 1
    # cds_annot['biotype'] = 'protein_coding'
    # cds_annot = cds_annot[['locus_tag', 'biotype', 'gene_name', 'gene_length']]
    # # export PC counts
    # quant_merged_pc = quant_merged[
    #     quant_merged["feature_id"].isin(cds_annot["locus_tag"].tolist())
    # ]
    quant_merged_pc = quant_merged
    outf2 = os.path.join(out_dir, "gene_counts_pc.tsv")
    quant_merged_pc.to_csv(outf2, index=False, sep="\t")
    # export ref gene df (using kallisto gene lengths instead of GFF)
    # cds_annot = pd.merge(cds_annot, gene_lengths, on="locus_tag")
    # cds_annot = cds_annot[["locus_tag", "biotype", "gene_name", "length"]]
    # cds_annot = cds_annot.drop_duplicates(keep="first")
    # cds_annot = cds_annot.rename(columns={"length": "gene_length"})
    # when not using GFF annotation:
    gene_lengths = gene_lengths.rename(columns={"length": "gene_length"})
    outf3 = os.path.join(out_dir, "ref_gene_df.tsv")
    # cds_annot.to_csv(outf3, index=False, sep="\t")
    gene_lengths.to_csv(outf3, index=False, sep="\t")


if __name__ == "__main__":
    parse()
