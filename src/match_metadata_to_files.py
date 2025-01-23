#!/usr/bin/env python
"""Read sample IDs and names from a TSV and match to filenamess in the given folder.

Print matched sample ID, name, and file path to stdout, separated by tabs.

Only emit file paths matching a given sample ID; if a file in the folder matches the
extension glob (.bam/fastq(.gz)), but does not match any sample ID in the input TSV,
then it will be skipped.
"""

import argparse
import csv
import sys
from pathlib import Path

AP = argparse.ArgumentParser(__doc__)
AP.add_argument("sample_in_metadata", type=Path)
AP.add_argument("sample_folder", type=Path)
args = AP.parse_args()


VALID_EXTENSIONS = (".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz")
available_files = [
    f for f in Path(args.sample_folder).iterdir() if f.name.endswith(VALID_EXTENSIONS)
]
print("Available files in", args.sample_folder, ":\n", available_files, file=sys.stderr)

with Path(args.sample_in_metadata).open() as infile:
    reader = csv.DictReader(infile, dialect="excel-tab")
    line_no = -1
    for line_no, row in enumerate(reader):  # noqa: B007
        sample_id = row.get("sample_unique_id")
        sample_name = row.get("sample_display_name")
        if not sample_id or not sample_name:
            sys.exit(
                "Error: The given metadata TSV file is missing required columns "
                "'sample_unique_id' or 'sample_display_name'."
            )

        matching_files = [f for f in available_files if sample_id in f.name]
        if len(matching_files) != 1:
            sys.exit(
                f"Error: sample_unique_id '{sample_id}' matches {len(matching_files)} "
                "files. Expected exactly one match."
            )

        seq_reads_file = matching_files[0]
        print(sample_id, sample_name, seq_reads_file.name, sep="\t")

    if line_no == -1:
        sys.exit("Error: The TSV file is empty or contains no valid data.")
    print(f"Loaded {line_no + 1} sample paths and identifiers", file=sys.stderr)
