#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sample_file",
        help = "tsv file containing sample metadata, with the following column " +
        "headers: sample, file1, file2, group, rep_no"
    )
    parser.add_argument(
        "data_dir", help = "Directory containing RNA-seq data"
    )
    parser.add_argument("outf", help="File for results")
    args = parser.parse_args()
    make_meta(**vars(args))

def make_meta(sample_file, data_dir, outf):
    sample_dat = pd.read_csv(sample_file, sep = "\t")
    # remove empty rows
    sample_dat.replace("", float("NaN"), inplace=True)
    sample_dat.dropna(subset = ["sample"], inplace=True)
    if(len(sample_dat.file2.value_counts()) > 0):
        f2_bnames = [os.path.basename(f) for f in sample_dat['file2'].to_list()]
        f2_full_paths = [os.path.join(data_dir, b) for b in f2_bnames]
        sample_dat['file2'] = f2_full_paths
        sample_dat['paired'] = "1"
    else:
        sample_dat['paired'] = "0"
    # remove empty columns 
    sample_dat.dropna(axis=1, how='all', inplace=True)
    # ensure that paths to files are correct
    f1_bnames = [os.path.basename(f) for f in sample_dat['file1'].to_list()]
    f1_full_paths = [os.path.join(data_dir, b) for b in f1_bnames]
    sample_dat['file1'] = f1_full_paths
    sample_dat.to_csv(outf, index=False, sep='\t')

if __name__ == "__main__":
    parse()
