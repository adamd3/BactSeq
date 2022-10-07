#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sample_file",
        help = "tsv file containing sample metadata, with the following column " +
        "headers: sample, file1, file2, group, rep_no, paired"
    )
    parser.add_argument(
        "data_dir", help = "Directory containing RNA-seq data"
    )
    parser.add_argument("outf", help="File for results")
    args = parser.parse_args()
    merge_meta(**vars(args))

def merge_meta(sample_file, data_dir, outf):
    sample_dat = pd.read_csv(sample_file, sep = "\t")
    # sample_dat['path_to_file1'] = [os.path.join(data_dir, s) for s in sample_dat['file1'].tolist()]
    # sample_dat['path_to_file2'] = [os.path.join(data_dir, s) for s in sample_dat['file2'].tolist()]
    sample_dat.to_csv(outf, index=False, sep='\t')

if __name__ == "__main__":
    parse()
