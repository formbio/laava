#!/usr/bin/env python3
"""Write sample ID, name, sequencing run ID, and optionally version as a 2-line TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pysam


def seq_run_id_from_bam(bam_fname):
    """Extract the sequencing run ID from the first read ID."""
    with pysam.AlignmentFile(bam_fname, "rb", check_sq=False) as bam_reader:
        first_read = next(iter(bam_reader))
        return first_read.qname.split("/")[0]


if __name__ == "__main__":
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("sample_id")
    AP.add_argument("sample_name")
    AP.add_argument("mapped_reads")
    AP.add_argument("-v", "--version", help="Workflow version string")
    AP.add_argument("-o", "--output", required=True, help="*.tsv")
    args = AP.parse_args()

    fieldnames = ["sample_unique_id", "sample_display_name", "sequencing_run_id"]
    outrow = [args.sample_id, args.sample_name, seq_run_id_from_bam(args.mapped_reads)]
    if args.version:
        fieldnames.append("version")
        outrow.append(args.version)

    with Path(args.output).open("w") as out_handle:
        cw = csv.writer(out_handle, delimiter="\t")
        cw.writerow(fieldnames)
        cw.writerow(outrow)
