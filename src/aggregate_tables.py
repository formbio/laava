#!/usr/bin/env python
"""Aggregate raw data tables into summary tables for the final report."""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def main(args):
    """CLI entry point to the analysis.

    Steps:

    1. Load the raw data table TSVs from the input path prefix
    2. Calculate aggregate summary tables from the loaded raw data.
    3. Output new TSVs (not .gz).  Optionally, also log to console.
    """
    # Run analysis
    results = analyze_alignments(args.path_prefix)
    # Write/show results
    write_results(results, args.path_prefix, args.verbose)


def analyze_alignments(path_prefix):
    """Analyze alignment data from TSV files.

    Parameters:
        path_prefix (str): Path prefix for input files

    Returns:
        dict: Dictionary containing analysis results and raw data
    """
    read_df = pd.read_csv(f"{path_prefix}.per_read.tsv.gz", sep="\t")
    # Filter to vector-only reads
    read_vector = read_df[read_df["reference_label"] == "vector"]

    total_read_count_all = read_df["effective_count"].sum()
    total_read_count_vector = read_vector["effective_count"].sum()

    result = {
        "agg_ref_type": get_ref_type_agg(
            read_df, total_read_count_all, total_read_count_vector
        ),
        "agg_subtype": get_subtype_agg(
            read_vector, total_read_count_all, total_read_count_vector
        ),
    }

    # Optional flip/flop aggregation
    ff_path = f"{path_prefix}.flipflop.tsv.gz"
    if os.path.exists(ff_path):
        ff_df = pd.read_csv(ff_path, sep="\t")
        result["agg_flipflop"] = get_flipflop_agg(ff_df)

    return result


def get_ref_type_agg(read_df, total_read_count_all, total_read_count_vector):
    """Counts and percentages of reference labels and types."""
    df = (
        read_df.groupby(["reference_label", "assigned_type"], dropna=False)[
            "effective_count"
        ]
        .sum()
        .reset_index(name="effective_count")
    )
    df = df.sort_values(
        ["reference_label", "effective_count"], ascending=[False, False]
    )
    pct_vector = round(df["effective_count"] * 100 / total_read_count_vector, 2)
    # Only show pct_vector for vector reads, otherwise NaN
    pct_vector[df["reference_label"] != "vector"] = np.nan
    df["pct_vector"] = pct_vector
    df["pct_total"] = round(df["effective_count"] * 100 / total_read_count_all, 2)
    return df


def get_subtype_agg(read_vector, total_read_count_all, total_read_count_vector):
    """Counts and percentages of assigned types and subtypes."""
    df = (
        read_vector.groupby(["assigned_type", "assigned_subtype"])["effective_count"]
        .sum()
        .reset_index(name="effective_count")
    )
    df = df.sort_values(["assigned_type", "effective_count"], ascending=[False, False])
    df["pct_vector"] = round(df["effective_count"] * 100 / total_read_count_vector, 2)
    df["pct_total"] = round(df["effective_count"] * 100 / total_read_count_all, 2)
    return df


def get_flipflop_agg(ff_df):
    """Counts of flip-flop configurations."""
    df = (
        ff_df.groupby(["type", "subtype", "leftITR", "rightITR"])
        .size()
        .reset_index(name="count")
    )
    return df


def write_results(results, path_prefix, verbose):
    """Save analysis results to CSV files."""
    # Save each analysis result
    for key, res_df in results.items():
        write_table(res_df, f"{path_prefix}.{key}.tsv", verbose)
    if verbose:
        print("\nResults saved to:", os.path.dirname(path_prefix))


def write_table(df, out_path, verbose):
    df.to_csv(out_path, index=False, sep="\t", float_format="%.2f")
    if verbose:
        print(out_path, ":")
        print(df.to_string(index=False, na_rep=""))
        print()
    else:
        print("Wrote", out_path, file=sys.stderr)


if __name__ == "__main__":
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument(
        "--path-prefix",
        type=str,
        required=True,
        help="""Path prefix for input files. The script will look for
        [prefix].alignments.tsv.gz and [prefix].per_read.tsv.gz""",
    )
    AP.add_argument(
        "--output-dir",
        type=str,
        default=".",
        help="Directory to save output files (default: current directory)",
    )
    AP.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print output TSV contents to stdout, too.",
    )
    main(AP.parse_args())
