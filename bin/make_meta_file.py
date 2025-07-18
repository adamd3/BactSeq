#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path
import sys


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sample_file",
        help="tsv file containing sample metadata, with the following column "
        "headers: sample, file1, file2, group, rep_no",
    )
    parser.add_argument("data_dir", help="Directory containing RNA-seq data")
    parser.add_argument("outf", help="File for results")
    args = parser.parse_args()
    make_meta(**vars(args))


def make_meta(sample_file, data_dir, outf):
    print(f"DEBUG: sample_file: {sample_file}", file=sys.stdout)
    print(f"DEBUG: data_dir: {data_dir}", file=sys.stdout)
    print(f"DEBUG: outf: {outf}", file=sys.stdout)

    if not os.path.exists(sample_file):
        print(f"ERROR: sample_file does not exist: {sample_file}", file=sys.stdout)
        sys.exit(1)

    try:
        sample_dat = pd.read_csv(sample_file, sep="\t")
    except Exception as e:
        print(f"ERROR: Could not read sample file: {e}", file=sys.stdout)
        sys.exit(1)

    print(f"DEBUG: Initial sample_dat shape: {sample_dat.shape}", file=sys.stdout)

    sample_dat = pd.read_csv(sample_file, sep="\t")

    # Remove empty rows
    sample_dat.replace("", float("NaN"), inplace=True)
    sample_dat.dropna(subset=["sample"], inplace=True)

    if len(sample_dat.file2.value_counts()) > 0:
        f2_bnames = [os.path.basename(f) for f in sample_dat["file2"].to_list()]
        f2_full_paths = [os.path.join(data_dir, b) for b in f2_bnames]
        sample_dat["file2"] = f2_full_paths
        sample_dat["paired"] = "1"
    else:
        sample_dat["paired"] = "0"

    # Remove empty columns
    sample_dat.dropna(axis=1, how="all", inplace=True)
    print(f"DEBUG: After dropna (columns): {sample_dat.shape}", file=sys.stdout)

    # Ensure that paths to files are correct
    f1_bnames = [os.path.basename(f) for f in sample_dat["file1"].to_list()]
    f1_full_paths = [os.path.join(data_dir, b) for b in f1_bnames]
    sample_dat["file1"] = f1_full_paths

    try:
        sample_dat.to_csv(outf, index=False, sep="\t")
        print(f"DEBUG: Successfully wrote to {outf}", file=sys.stdout)
        print(f"DEBUG: File {outf} exists: {os.path.exists(outf)}", file=sys.stdout)
        print(f"DEBUG: File {outf} size: {os.path.getsize(outf) if os.path.exists(outf) else 'N/A'}", file=sys.stdout)
    except Exception as e:
        print(f"ERROR: Could not write output file {outf}: {e}", file=sys.stdout)
        sys.exit(1)


if __name__ == "__main__":
    parse()
