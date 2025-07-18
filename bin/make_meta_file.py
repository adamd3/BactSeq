#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path
import sys
from urllib.parse import urljoin


def join_path_or_url(base, filename):
    """Join a path or URL with a filename appropriately."""
    if base.startswith(('http://', 'https://')):
        # For URLs, ensure base ends with / and use urljoin
        if not base.endswith('/'):
            base += '/'
        return urljoin(base, filename)
    else:
        # For local paths, use os.path.join
        return os.path.join(base, filename)


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
    if not os.path.exists(sample_file):
        print(f"ERROR: sample_file does not exist: {sample_file}", file=sys.stdout)
        sys.exit(1)

    try:
        sample_dat = pd.read_csv(sample_file, sep="\t")
    except Exception as e:
        print(f"ERROR: Could not read sample file: {e}", file=sys.stdout)
        sys.exit(1)

    sample_dat = pd.read_csv(sample_file, sep="\t")

    # Remove empty rows
    sample_dat.replace("", float("NaN"), inplace=True)
    sample_dat.dropna(subset=["sample"], inplace=True)

    if len(sample_dat.file2.value_counts()) > 0:
        f2_bnames = [os.path.basename(f) for f in sample_dat["file2"].to_list()]
        f2_full_paths = [join_path_or_url(data_dir, b) for b in f2_bnames]
        sample_dat["file2"] = f2_full_paths
        sample_dat["paired"] = "1"
    else:
        sample_dat["paired"] = "0"

    # Remove empty columns
    sample_dat.dropna(axis=1, how="all", inplace=True)

    # Ensure that paths to files are correct
    f1_bnames = [os.path.basename(f) for f in sample_dat["file1"].to_list()]
    f1_full_paths = [join_path_or_url(data_dir, b) for b in f1_bnames]
    sample_dat["file1"] = f1_full_paths

    try:
        sample_dat.to_csv(outf, index=False, sep="\t")
    except Exception as e:
        print(f"ERROR: Could not write output file {outf}: {e}", file=sys.stdout)
        sys.exit(1)


if __name__ == "__main__":
    parse()
