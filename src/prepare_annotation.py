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
from collections.abc import Iterable
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


def write_annotation_txt(out_fname: str, bed_rows: dict, other_rows: Iterable):
    """Write PacBio-style annotations to `out_fname`.

    Take the vector annotations and non-'vector' sequence names and labels, format
    it for annotation.txt, and append it to the same output file.

    Skip any duplicate 'vector' sequence label appearing in the other references, and
    catch if a sequence name is reused across multiple labels.
    """
    vector_row = bed_rows["vector"]
    repcap_row = bed_rows["repcap"]
    with Path(out_fname).open("w+") as outf:
        outf.write("NAME={};TYPE={};REGION={}-{};\n".format(*vector_row))
        seen_seq_names_and_labels = {vector_row.seq_name: vector_row.ref_label}
        if repcap_row:
            outf.write("NAME={};TYPE={};REGION={}-{};\n".format(*repcap_row))
            seen_seq_names_and_labels[repcap_row.seq_name] = repcap_row.ref_label
        for orow in other_rows:
            if orow.ref_label == "vector":
                if orow.seq_name != vector_row.seq_name:
                    raise RuntimeError(
                        "Reference label 'vector' listed for additional "
                        f"reference names, but sequence name {orow.seq_name} does not "
                        f"match the previously given {vector_row.seq_name}"
                    )
                continue
            if orow.seq_name in seen_seq_names_and_labels:
                prev_type = seen_seq_names_and_labels[orow.seq_name]
                if orow.ref_label != prev_type:
                    raise RuntimeError(
                        f"Sequence name {orow.seq_name} listed with "
                        f"different labels: first {prev_type}, then "
                        f"{orow.ref_label}"
                    )
                continue
            outf.write("NAME={seq_name};TYPE={ref_label};\n".format(**orow._asdict()))


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
    AP.add_argument("-o", "--output", help="Output filename (*.txt).")
    args = AP.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    try:
        bed_rows = read_annotation_bed(args.annotation_bed, args.itr_labels)
        otr_rows = read_reference_names(args.reference_names)
        write_annotation_txt(args.output, bed_rows, otr_rows)
    except RuntimeError as exc:
        sys.exit(str(exc))
