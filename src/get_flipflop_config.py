#!/usr/bin/env python3
"""Get ITR flip flop configurations.

Must have already run `summarize_AAV_alignment.py` to get a .tagged.BAM file!
"""

import sys
from csv import DictReader
from typing import NamedTuple

import parasail
import pysam
from Bio import SeqIO


SW_SCORE_MATRIX = parasail.matrix_create("ACGT", 2, -5)

SEQ_AAV2 = dict(
    left_flip="TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT",
    left_flop="TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT",
    right_flip="AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA",
    right_flop="AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA",
)


class FlipFlopSeqSet(NamedTuple):
    """Immutable container for flip-flop sequences."""

    left_flip: str
    left_flop: str
    right_flip: str
    right_flop: str

    @classmethod
    def from_fasta(cls, fname):
        flipflip_seqs = {
            "left_flip": None,
            "left_flop": None,
            "right_flip": None,
            "right_flop": None,
        }
        known_seq_ids = {
            "SEQ_LEFT_FLIP": "left_flip",
            "SEQ_LEFT_FLOP": "left_flop",
            "SEQ_RIGHT_FLIP": "right_flip",
            "SEQ_RIGHT_FLOP": "right_flop",
        }
        print(f"Reading {fname}....")
        for r in SeqIO.parse(fname, "fasta"):
            if r.id in known_seq_ids:
                flipflip_seqs[known_seq_ids[r.id]] = str(r.seq)
            else:
                print(
                    "WARNING: Sequence IDs must be "
                    "SEQ_LEFT_FLIP|SEQ_LEFT_FLOP|SEQ_RIGHT_FLIP|SEQ_RIGHT_FLOP. "
                    f"Is {r.id} instead. Ignoring!"
                )
        if not all(flipflip_seqs.keys()):
            missing_keys = [
                key for key, val in known_seq_ids.items() if flipflip_seqs[val] is None
            ]
            raise RuntimeError(
                f"Required flip-flop sequence IDs missing from {fname}: " ", ".join(
                    missing_keys
                )
            )
        return cls(**flipflip_seqs)


def identify_flip_flop(r, ff_seq):
    """Determine left and right flip/flip/unclassified configurations.

    Assume record tag:AT is vector, tag:AX can be full|left-partial|right-partial|partial
    Add back a tag 'AF' that is [flip/flop]-[flip/flop]
    """
    min_score = 250

    t = dict(r.tags)
    if t["AX"] not in (
        "vector-full",
        "vector-left-partial",
        "vector-right-partial",
        "vector-partial",
    ):
        raise RuntimeError(
            "Input BAM records must have a `AX` tag assigned by first running "
            "summarize_AAV_alignment.py. Abort!"
        )

    config_left, config_right = "unclassified", "unclassified"
    if t["AX"] == "vector-partial":
        # ignore, since both sides are missing chunks of ITR
        return "unclassified", "unclassified"

    if t["AX"] in ("vector-full", "vector-left-partial"):
        o1 = parasail.sw_trace(r.query, ff_seq.left_flip, 3, 1, SW_SCORE_MATRIX)
        o2 = parasail.sw_trace(r.query, ff_seq.left_flop, 3, 1, SW_SCORE_MATRIX)
        if o1.score > o2.score and o1.score > min_score:
            config_left = "flip"
        elif o2.score > o1.score and o2.score > min_score:
            config_left = "flop"
        else:
            config_left = "unclassified"

    if t["AX"] in ("vector-full", "vector-right-partial"):
        o1 = parasail.sw_trace(
            r.query[-len(ff_seq.right_flip) - 10 :],
            ff_seq.right_flip,
            3,
            1,
            SW_SCORE_MATRIX,
        )
        o2 = parasail.sw_trace(
            r.query[-len(ff_seq.right_flop) - 10 :],
            ff_seq.right_flop,
            3,
            1,
            SW_SCORE_MATRIX,
        )
        if o1.score > o2.score and o1.score > min_score:
            config_right = "flip"
        elif o2.score > o1.score and o2.score > min_score:
            config_right = "flop"
        else:
            config_right = "unclassified"
    return config_left, config_right


def load_per_read_info(fname):
    """Load per-read info, keyed by read IDs, from a CSV file."""
    with open(fname) as in_csv:
        read_info = {r["read_id"]: r for r in DictReader(in_csv, delimiter="\t")}
    return read_info


def main(per_read_csv, tagged_bam, output_prefix, flipflop_fasta):
    """Entry point."""
    if flipflop_fasta is None:
        flipflop_seqs = FlipFlopSeqSet(**SEQ_AAV2)
    else:
        flipflop_seqs = FlipFlopSeqSet.from_fasta(flipflop_fasta)

    read_info = load_per_read_info(per_read_csv)

    with open(output_prefix + ".flipflop_assignments.txt", "w") as fout:
        fout.write("name\ttype\tsubtype\tstart\tend\tleftITR\trightITR\n")
        reader = pysam.AlignmentFile(open(tagged_bam), "rb", check_sq=False)
        writer1 = pysam.AlignmentFile(
            open(output_prefix + ".vector-full-flipflop.bam", "w"),
            "wb",
            header=reader.header,
        )
        writer2 = pysam.AlignmentFile(
            open(output_prefix + ".vector-leftpartial-flipflop.bam", "w"),
            "wb",
            header=reader.header,
        )
        writer3 = pysam.AlignmentFile(
            open(output_prefix + ".vector-rightpartial-flipflop.bam", "w"),
            "wb",
            header=reader.header,
        )
        for r in reader:
            t = dict(r.tags)
            if t["AT"] == "vector" and t["AX"] in (
                "vector-full",
                "vector-left-partial",
                "vector-right-partial",
            ):
                c_l, c_r = identify_flip_flop(r, flipflop_seqs)
                d = r.to_dict()
                a_type = read_info[r.qname]["assigned_type"]
                if a_type not in ("scAAV", "ssAAV"):
                    continue
                d["tags"].append("AF:Z:" + c_l + "-" + c_r)
                d["tags"].append("AG:Z:" + a_type)
                if t["AX"] == "vector-full":
                    writer = writer1
                elif t["AX"] == "vector-right-partial":
                    writer = writer2
                elif t["AX"] == "vector-left-partial":
                    writer = writer3
                writer.write(pysam.AlignedSegment.from_dict(d, r.header))
                fout.write(
                    "\t".join(
                        [
                            r.qname,
                            a_type,
                            t["AX"],
                            str(r.reference_start),
                            str(r.reference_end),
                            c_l,
                            c_r,
                        ]
                    )
                    + "\n"
                )

        writer1.close()
        writer2.close()
        writer3.close()
        print("Output summmary:", fout.name)

    print(
        f"Individual BAM files written: {output_prefix}.vector- full,leftpartial,rightpartial -flipflop.bam"
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("sorted_tagged_bam", help="Sorted tagged BAM file")
    parser.add_argument("per_read_csv", help="Per read CSV file")
    parser.add_argument("-o", "--output-prefix", help="Output prefix", required=True)
    parser.add_argument(
        "--flipflop-fasta",
        default=None,
        help="(optional) flip flop fasta file (if not given, uses AAV2 default)",
    )

    args = parser.parse_args()

    main(
        args.per_read_csv,
        args.sorted_tagged_bam,
        args.output_prefix,
        args.flipflop_fasta,
    )
