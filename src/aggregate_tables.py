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
    alignments_df = pd.read_csv(f"{path_prefix}.alignments.tsv.gz", sep="\t")
    
    # Pre-aggregate the alignment data by read_id (sum map_len)
    map_len_df = alignments_df[alignments_df["is_mapped"] == "Y"].groupby("read_id")["map_len"].sum().reset_index()
    
    # Filter to vector-only reads
    read_vector = read_df[read_df["reference_label"] == "vector"]

    # Filter to lambda-only reads
    read_lambda = read_df[(read_df["reference_label"] == "Lambda") | (read_df["reference_label"] == "lambda")]

    total_read_count_all = read_df["effective_count"].sum()
    total_read_count_vector = read_vector["effective_count"].sum()
    total_read_count_lambda = 0 if read_lambda.empty else read_lambda["effective_count"].sum()

    result = {
        "agg_ref_type": get_ref_type_agg(
            read_df, map_len_df, total_read_count_all, total_read_count_vector, total_read_count_lambda
        ),
        "agg_subtype": get_subtype_agg(
            read_vector, map_len_df, total_read_count_all, total_read_count_vector, total_read_count_lambda
        ),
    }

    # Optional flip/flop aggregation
    ff_path = f"{path_prefix}.flipflop.tsv.gz"
    if os.path.exists(ff_path):
        ff_df = pd.read_csv(ff_path, sep="\t")
        result["agg_flipflop"] = get_flipflop_agg(ff_df)

    return result


def get_ref_type_agg(read_df, map_len_df, total_read_count_all, total_read_count_vector, total_read_count_lambda):
    """Counts and percentages of reference labels and types."""
    # Merge read_df with pre-aggregated map_len data
    merged_df = pd.merge(
        read_df[["read_id", "reference_label", "assigned_type", "effective_count"]], 
        map_len_df,
        on="read_id",
        how="left"  # Keep all reads, even if they don't have alignments
    )
    merged_df["map_len"] = merged_df["map_len"].fillna(0)
    
    # Do a single groupby to get both effective_count and base
    df = merged_df.groupby(["reference_label", "assigned_type"], dropna=False).agg({
        "effective_count": "sum",
        "map_len": "sum"
    }).reset_index()
    
    # Rename map_len to base
    df = df.rename(columns={"map_len": "base"})
    
    df = df.sort_values(
        ["reference_label", "effective_count"], ascending=[False, False]
    )
    pct_vector = round(df["effective_count"] * 100 / total_read_count_vector, 2)
    # Only show pct_vector for vector reads, otherwise NaN
    pct_vector[df["reference_label"] != "vector"] = np.nan
    df["pct_vector"] = pct_vector
    df["pct_total"] = round(df["effective_count"] * 100 / total_read_count_all, 2)
    df["counts_wo_lambda"] = df["effective_count"]
    df["pct_wo_lambda"] = round(df["effective_count"] * 100 / (total_read_count_all-total_read_count_lambda), 2)
    df.loc[df["reference_label"].isin(["lambda", "Lambda"]), "counts_wo_lambda"] = 0
    df.loc[df["reference_label"].isin(["lambda", "Lambda"]), "pct_wo_lambda"] = 0.0
    return df


def get_subtype_agg(read_vector, map_len_df, total_read_count_all, total_read_count_vector, total_read_count_lambda):
    """Counts and percentages of assigned types and subtypes."""
    # Handle case where there are no vector reads
    if read_vector.empty:
        # Return DataFrame with NA values that R can handle without affecting calculations
        # Using pd.NA ensures R recognizes this as missing data that won't interfere with sums
        return pd.DataFrame({
            "assigned_type": [pd.NA],
            "assigned_subtype": [pd.NA],
            "effective_count": [0],  # Use 0 so R sum() works correctly
            "base": [0],
            "pct_vector": [0.0],
            "pct_total": [0.0], 
            "pct_wo_lambda": [0.0]
        })
    
    # Merge read_vector with pre-aggregated map_len data
    merged_df = pd.merge(
        read_vector[["read_id", "assigned_type", "assigned_subtype", "effective_count"]],
        map_len_df,
        on="read_id",
        how="left"
    )
    merged_df["map_len"] = merged_df["map_len"].fillna(0)
    
    # Group by to get both effective_count and base
    df = merged_df.groupby(["assigned_type", "assigned_subtype"]).agg({
        "effective_count": "sum",
        "map_len": "sum"
    }).reset_index()
    
    # Rename map_len to base
    df = df.rename(columns={"map_len": "base"})
    
    df = df.sort_values(["assigned_type", "effective_count"], ascending=[False, False])
    df["pct_vector"] = round(df["effective_count"] * 100 / total_read_count_vector, 2)
    df["pct_total"] = round(df["effective_count"] * 100 / total_read_count_all, 2)
    df["pct_wo_lambda"] = round(df["effective_count"] * 100 / (total_read_count_all-total_read_count_lambda), 2)
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
