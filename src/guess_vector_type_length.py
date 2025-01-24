#!/usr/bin/env python3
"""Infer single stranded (ssAAV) vs. self-complementary (scAAV) vector type.

Approach:
- scAAV packages best when the expression cassette (ITR to ITR) is 2.3-2.4 kb long.
- ssAAV constructs are prone to unresolved dimers when the length is much less than
  4.7kb; the drug product is generally poor for lengths below 4.0kb and stuffer
  sequences are necessary to rescue the construct.
- Therefore, we read the length of the vector's annotated ITR-to-ITR region to determine
  the length of the expression cassette. If this length is below a threshold, we report
  scAAV, otherwise ssAAV.
- If the heuristic guesses incorrectly, the mistake will be prominent in the report:
  most vector reads will be marked as "other-vector" instead of "scAAV" or "ssAAV".
"""

from __future__ import annotations

import argparse
import sys

from prepare_annotation import read_annotation_bed

SC_MAX_THRESHOLD = 2300


def length_from_annotation(fname: str, itr_labels: list[str]) -> int:
    """Read annotation BED to calculate the expression cassette length."""
    vec = read_annotation_bed(fname, itr_labels)["vector"]
    return vec.end - vec.start1 + 1


if __name__ == "__main__":
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("ann_file", help="Vector annotation (.bed)")
    AP.add_argument("itr_labels", nargs="*", help="ITR label(s) in annotation BED")
    AP.add_argument(
        "-t",
        "--sc-max-threshold",
        default=SC_MAX_THRESHOLD,
        type=int,
        help="""Maximum length of the expression cassette, including ITRs, in the
        given vector annotation to report as scAAV. [Default: %(default)d]""",
    )
    args = AP.parse_args()
    cassette_len = length_from_annotation(args.ann_file, args.itr_labels)
    print("Expression cassette length is", cassette_len, file=sys.stderr)
    if cassette_len <= args.sc_max_threshold:
        print("sc")
    else:
        print("ss")
