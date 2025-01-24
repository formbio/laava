#!/usr/bin/env python3
"""Format vector annotations and other references as PacBio annotation.txt.

The output annotation.txt is a bespoke format required by the PacBio scripts. It looks
like::

    NAME=myVector;TYPE=vector;REGION=1795-6553;
    NAME=myCapRep;TYPE=repcap;
    NAME=myHelper;TYPE=helper;
    NAME=chr1;TYPE=host;
    NAME=chr2;TYPE=host;

Or, in a dual construct (repcap in vector backbone)::

    NAME=myVector;TYPE=vector;REGION=1795-6553;
    NAME=myVector;TYPE=repcap;REGION=1895-5987;
    NAME=myHelper;TYPE=helper;
    NAME=chr1;TYPE=host;
    NAME=chr2;TYPE=host;


Specifically:

  - NAME is a sequence name/ID used in the reference sequence files that reads will be
    mapped to.
  - TYPE is the source "type" or label of the reference sequence -- one of 'vector', 'helper',
    'repcap', 'host', or 'lambda'.
  - There must be exactly one line where TYPE is 'vector'.
  - The 'vector' line also has a REGION field with the 1-based start and end positions.
    In AAV, this is the payload region including the ITRs.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import NamedTuple


class AnnRow(NamedTuple):
    seq_name: str
    ref_label: str
    start1: int
    end: int


# From summarize_alignment.py
KNOWN_ANNOT_LABELS = {"vector", "repcap", "helper", "lambda", "host"}


def read_annotation_bed(fname: str, itr_labels: list[str]):
    """Read a UCSC BED4 file and extract the coordinates of the ITR-to-ITR region.

    Scan a UCSC BED4 file for the region coordinates of (a) the inverted terminal
    repeats (ITRs), identified by `itr_labels`, or if those are not provided, (b) a row
    with the label exactly named 'vector' (legacy mode).

    Check that the input contains exactly (a) two ITR rows, or (b) one 'vector' row.
    Ignore all other BED rows.
    """
    out_rows = {
        "repcap": None,
        "vector": None,
    }
    if len(itr_labels):
        itr_labels = set(filter(None, itr_labels))
    itr_slots = []

    # Collect known annotations from the BED file
    with Path(fname).open() as infile:
        for line in infile:
            # Require BED4 or more
            seq_name, start0, end, label = line.rstrip().split("\t")[:4]
            region = AnnRow(seq_name, label, int(start0) + 1, int(end))
            if label in ("vector", "repcap"):
                # Only allow 1 of each
                if out_rows[label] is not None:
                    raise RuntimeError(
                        f"Input {args.annotation_bed} contains more than one row "
                        f"labeled '{label}'."
                    )
                out_rows[label] = region
            elif label in itr_labels:
                # Assign to left and then right slot; labels 1 & 2 are interchangeable
                if len(itr_slots) > 1:
                    raise RuntimeError(
                        f"Input {args.annotation_bed} contains more than two rows "
                        f"with potential ITR labels ('{itr_labels}')."
                    )
                if len(itr_slots) == 1 and itr_slots[0].seq_name != seq_name:
                    raise RuntimeError(
                        f"First observed ITR label {itr_slots[0].label} "
                        f"is on sequence '{itr_slots[0].seq_name}', but second ITR "
                        f"label {label} is on a different sequence '{seq_name}'."
                    )
                itr_slots.append(region)
                logging.info("Found an ITR region: %s", region)
            else:
                logging.info("Ignoring unrecognized label: %s", label)

    # Verify the recognized annotations
    if itr_labels:
        orig_vector = out_rows.get("vector")
        if len(itr_slots) == 2:
            new_vector = AnnRow(
                itr_slots[0].seq_name,
                "vector",
                min(itr_slots[0].start1, itr_slots[1].start1),
                max(itr_slots[0].end, itr_slots[1].end),
            )
            out_rows["vector"] = new_vector
            # Use this as "vector"
            if orig_vector:
                logging.info(
                    "Taking vector coordinates (%s, %s, %s) from ITR labels, "
                    "overriding explicit 'vector' coordinates (%s, %s, %s)",
                    new_vector.seq_name,
                    new_vector.start1,
                    new_vector.end,
                    orig_vector.seq_name,
                    orig_vector.start1,
                    orig_vector.end,
                )
        elif len(itr_slots) == 1:
            raise RuntimeError(
                f"ITR labels were specified as {itr_labels}; expected to find 2, "
                f"but only found 1 ({itr_slots[0]}) in "
                f"the input annotation file {args.annotation_bed}"
            )
        else:
            raise RuntimeError(
                f"ITR labels were specified as {itr_labels}, but were not found in "
                f"the input annotation file {args.annotation_bed}"
            )
    elif out_rows["vector"] is None:
        # Legacy mode: require 'vector' label if ITR labels were not specified
        raise RuntimeError(
            f"Input {args.annotation_bed} must contain a row labeled 'vector'."
        )

    return out_rows


def read_reference_names(fname: str):
    """Read a 2-column TSV of reference sequence names and labels."""
    with Path(fname).open() as infile:
        # Skip header
        lines = iter(infile)
        next(infile)
        for line in lines:
            seq_name, ref_label = line.split()
            if ref_label not in KNOWN_ANNOT_LABELS:
                logging.info(
                    "Nonstandard reference label %s; "
                    "the known labels are: vector, repcap, helper, host",
                    ref_label,
                )
            yield AnnRow(seq_name, ref_label, None, None)


if __name__ == "__main__":
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument(
        "annotation_bed",
        help="Vector sequence annotations, as UCSC BED4 with 0-based coordinates.",
    )
    AP.add_argument(
        "reference_names",
        help="Reference sequence names and their labels, in 2 columns.",
    )
    AP.add_argument("itr_labels", nargs="*", help="ITR label(s) in annotation BED")
    args = AP.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    try:
        # Just exercise the parsing code
        bed_rows = read_annotation_bed(args.annotation_bed, args.itr_labels)
        otr_rows = read_reference_names(args.reference_names)
    except RuntimeError as exc:
        sys.exit(str(exc))
