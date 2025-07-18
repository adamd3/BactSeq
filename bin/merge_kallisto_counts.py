#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_f", help="TSV file containing metadata")
    parser.add_argument("--out_dir", help="Directory for results")
    args = parser.parse_args()
    merge_counts(**vars(args))


def merge_counts(metadata_f, out_dir):
    # Read metadata and merge counts
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
    
    # Get gene lengths from kallisto results
    quant_merged = quant_merged.rename_axis("feature_id").reset_index()
    gene_lengths = quant_merged[["feature_id", "length"]]
    gene_lengths = gene_lengths.rename(columns={"feature_id": "locus_tag"})
    quant_merged = quant_merged[["feature_id"] + metadata["sample"].tolist()]
    
    # Export merged counts
    outf1 = os.path.join(out_dir, "gene_counts.tsv")
    quant_merged.to_csv(outf1, index=False, sep="\t")
    
    # Export protein-coding counts (currently same as all counts)
    quant_merged_pc = quant_merged
    outf2 = os.path.join(out_dir, "gene_counts_pc.tsv")
    quant_merged_pc.to_csv(outf2, index=False, sep="\t")
    
    # Export reference gene dataframe with gene lengths
    gene_lengths = gene_lengths.rename(columns={"length": "gene_length"})
    outf3 = os.path.join(out_dir, "ref_gene_df.tsv")
    gene_lengths.to_csv(outf3, index=False, sep="\t")


if __name__ == "__main__":
    parse()
