#!/usr/bin/env python
"""Aggregate raw data tables into summary tables for the final report.

Details:
    1. Load the raw data table TSVs from the input path prefix
    2. Calculate aggregate summary tables from the loaded raw data.
    3. Output new TSVs (not .gz).  Optionally, also log to console.

Todo here:
    - rename to "aggregate_tables.py"
    - factor out the calculation of each table into a separate function (as possible)
    - ensure the assigned-types table shows all rows, not just vector
    - fix output TSV names
    - check that all summary tables of interest in the report are handled here

Separately:
    - reorganize the outputs of both Nextflow processes into folders
    - orchestrate this script in Nextflow make_report
        - insert as a step before create_report.R
    - update report.Rmd to load these new TSVs instead of computing them

"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze alignment data from TSV files."
    )
    parser.add_argument(
        "--path-prefix",
        type=str,
        required=True,
        help="""Path prefix for input files. The script will look for
        [prefix].alignments.tsv.gz and [prefix].per_read.tsv.gz""",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=".",
        help="Directory to save output files (default: current directory)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print output TSV contents to stdout, too.",
    )
    return parser.parse_args()


def analyze_alignments(path_prefix):
    """Analyze alignment data from TSV files.

    Parameters:
        path_prefix (str): Path prefix for input files

    Returns:
        dict: Dictionary containing analysis results and raw data
    """
    read_df = pd.read_csv(
        f"{path_prefix}.per_read.tsv.gz", sep="\t"
    )  # prev. x_all_read
    total_read_count_all = read_df["effective_count"].sum()

    # Filter to vector-only reads
    read_vector = read_df[
        read_df["reference_label"] == "vector"
    ].copy()  # was: x_read_vector
    total_read_count_vector = read_vector["effective_count"].sum()

    ref_type_agg = get_ref_type_agg(
        read_df, total_read_count_all, total_read_count_vector
    )
    subtype_agg = get_subtype_agg(
        read_vector, total_read_count_all, total_read_count_vector
    )

    return {
        "ref_type_agg": ref_type_agg,
        "subtype_agg": subtype_agg,
    }


def get_ref_type_agg(read_df, total_read_count_all, total_read_count_vector):
    """Get overall read distribution by reference and type."""
    ref_type_agg = (
        read_df
        .groupby(["reference_label", "assigned_type"], dropna=False)
        ["effective_count"].sum()
        .reset_index(name="e_count")
    )
    ref_type_agg = ref_type_agg.sort_values(
        ["reference_label", "e_count"], ascending=[False, False]
    )
    ref_type_agg["pct_total"] = round(
        ref_type_agg["e_count"] * 100 / total_read_count_all, 2
    )
    pct_vector = round(ref_type_agg["e_count"] * 100 / total_read_count_vector, 2)
    # Only show pct_vector for vector reads, otherwise NaN
    pct_vector[ref_type_agg["reference_label"] != "vector"] = np.nan
    ref_type_agg["pct_vector"] = pct_vector
    # was: df_read1, df_read_vector1
    return ref_type_agg


def get_subtype_agg(read_vector, total_read_count_all, total_read_count_vector):
    """Type and subtype analysis"""
    subtype_agg = (
        read_vector.groupby(["assigned_type", "assigned_subtype"])["effective_count"]
        .sum()
        .reset_index(name="e_count")
    )
    subtype_agg = subtype_agg.sort_values(
        ["assigned_type", "e_count"], ascending=[False, False]
    )
    subtype_agg["pct_vector"] = round(
        subtype_agg["e_count"] * 100 / total_read_count_vector, 2
    )
    subtype_agg["pct_total"] = round(
        subtype_agg["e_count"] * 100 / total_read_count_all, 2
    )
    # was: df_read_vector2, sopt3
    return subtype_agg


def write_results(results, output_dir, sample_id, verbose):
    """Save analysis results to CSV files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Save each analysis result
    write_table(
        results["ref_type_agg"], output_dir / f"{sample_id}.ref_type_agg.tsv", verbose
    )
    write_table(
        results["subtype_agg"], output_dir / f"{sample_id}.subtype_agg.tsv", verbose
    )


def write_table(df, out_path, verbose):
    df.to_csv(out_path, index=False)
    if verbose:
        print(out_path, ":")
        print(df.to_string(index=False, na_rep=""))
        print()
    else:
        print("Wrote", out_path, file=sys.stderr)


def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_args()

    # Run analysis
    results = analyze_alignments(args.path_prefix)
    sample_id = os.path.basename(args.path_prefix)

    # Write/show results
    write_results(results, args.output_dir, sample_id, args.verbose)
    print(f"\nResults saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
