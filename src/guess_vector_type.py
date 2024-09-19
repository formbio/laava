#!/usr/bin/env python3
"""Infer single stranded (ssAAV) vs. self-complementary (scAAV) vector type.

Approach:
- scAAV reads aligned to the vector contig almost always have two alignments (forward
  and reverse-complement), with rare single alignments coming from contaminants and the
  plasmid backbone.
- ssAAV reads aligned to the vector more often have single alignments; reads with two
  alignments most often come from "snapback" vector genomes.
- Therefore, we calculate the fraction of reads aligned to the vector contig that have
  exactly one alignment. If this fraction is below a threshold, we report scAAV,
  otherwise ssAAV.
"""

from __future__ import annotations

import argparse
import collections
import sys

import pysam

SC_MAX_THRESHOLD = 0.1


def count_alignments(sam_fname, vector_name):
    """Collect each read's length and aligned suffix counts."""

    def filter_records():
        reader = pysam.AlignmentFile(sam_fname, "r", check_sq=False)
        for rec in reader:
            if not rec.is_mapped:
                continue
            if vector_name and rec.reference_name != vector_name:
                continue
            _movie_id, read_num, _suffix = rec.qname.split("/", 2)
            yield read_num

    n_aln_per_read = collections.Counter(filter_records())
    histo_n_aln = collections.Counter(n_aln_per_read.values())
    print(histo_n_aln, file=sys.stderr)
    n_single = histo_n_aln[1]
    total = sum(histo_n_aln.values())
    frac_single = n_single / total
    print(f"{n_single} / {total} = {frac_single} singly-aligned reads", file=sys.stderr)
    return frac_single


if __name__ == "__main__":
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("sam_file", help="Read alignment in SAM format.")
    AP.add_argument("-v", "--vector-name", help="Vector reference sequence name.")
    AP.add_argument(
        "-t",
        "--sc-max-threshold",
        default=SC_MAX_THRESHOLD,
        type=float,
        help="""Maximum fraction of singly-aligned reads for this sample to be reported
        as self-complementary AAV. [Default: %(default)f]""",
    )
    args = AP.parse_args()
    frac = count_alignments(args.sam_file, args.vector_name)
    if frac <= args.sc_max_threshold:
        print("sc")
    else:
        print("ss")
