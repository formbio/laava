#!/usr/bin/env python3
"""Map reference sequence IDs to "reference label" categories for reporting.

Output is a 2-column TSV with no header:
1. FASTA sequence ID
2. Reference label, one of: "vector", "repcap", "helper", "host", "lambda", or the original
   sequence ID if it's not one of those.
"""

import argparse
import logging
from pathlib import Path

from Bio import SeqIO


AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument("vector", help="Vector plasmid reference sequence (FASTA).")
AP.add_argument(
    "--packaging", help="Packaging plasmid and other reference sequences (FASTA)."
)
AP.add_argument("--host", help="Host cell genome reference sequence (FASTA).")
AP.add_argument("--repcap-name", help="RepCap sequence ID in the packaging file.")
AP.add_argument("--helper-name", help="Helper sequence ID in the packaging file.")
AP.add_argument("--lambda-name", help="Lambda sequence ID in the packaging file.")
AP.add_argument("-o", "--output", required=True, help="Output filename.")


def _main(args):
    out_rows = []

    # Vector sequence -- should be just 1 in the FASTA; label it "vector"
    vector_rec = SeqIO.read(args.vector, "fasta")
    out_rows.append((vector_rec.name, "vector"))

    # Packaging sequences -- multi-FASTA; apply renaming
    if args.packaging:
        packaging_name_map = {
            seq_id: reference_label
            for seq_id, reference_label in [
                (args.repcap_name, "repcap"),
                (args.helper_name, "helper"),
                (args.lambda_name, "lambda"),
            ]
            if seq_id is not None
        }
        for pkg_rec in SeqIO.parse(args.packaging, "fasta"):
            label = packaging_name_map.get(pkg_rec.name, pkg_rec.name)
            out_rows.append((pkg_rec.name, label))

    # Host sequences -- multi-FASTA; label them all "host"
    if args.host:
        for host_rec in SeqIO.parse(args.host, "fasta"):
            out_rows.append((host_rec.name, "host"))

    # Write to file
    with Path(args.output).open("w") as outf:
        for col1, col2 in out_rows:
            outf.write(f"{col1}\t{col2}\n")
    logging.info("Wrote %s", args.output)


if __name__ == "__main__":
    args = AP.parse_args()
    _main(args)
