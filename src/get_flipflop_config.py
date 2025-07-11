#!/usr/bin/env python3
"""Get ITR flip flop configurations.

Must have already run `summarize_alignment.py` to get a .tagged.BAM file!
"""

from __future__ import annotations

import csv
import gzip
from typing import NamedTuple

import parasail
import pysam
from Bio import SeqIO
import re


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


def identify_flip_flop(df, ff_seq, vector_type="ss", orientation="left"):
    """Determine left and right flip/flip/unclassified configurations.

    Assume record tag:AT is vector, tag:AX can be full|left-partial|right-partial|partial
    Add back a tag 'AF' that is [flip/flop]-[flip/flop]
    """
    min_score = 250
    min_insert = 0 
    config_left, config_right = "unclassified", "unclassified"

    # Expected flags
    # 0 read forward strand
    # 16 read reverse strand
    # 2048 read forward strand + supplementary alignment
    # 2064 read reverse strand + supplementary alignment
    # 4 read unmapped

    if (vector_type == "sc"):
        if len(df) == 2:
            if orientation == "left":    
                forward_query = df[(df['flag'] == 0 ) | (df['flag'] == 2048)]['seq'].iloc[0]
                reverse_query = df[(df['flag'] == 16) | (df['flag'] == 2064)]['seq'].iloc[0]
                
                o11 = parasail.sw_trace(forward_query, ff_seq.left_flip, 3, 1, SW_SCORE_MATRIX)
                o12 = parasail.sw_trace(forward_query, ff_seq.left_flop, 3, 1, SW_SCORE_MATRIX)
                o21 = parasail.sw_trace(reverse_query, ff_seq.left_flip, 3, 1, SW_SCORE_MATRIX)
                o22 = parasail.sw_trace(reverse_query, ff_seq.left_flop, 3, 1, SW_SCORE_MATRIX)

                if o11.score > o12.score and o11.score > min_score:
                    config_left = "flip"
                elif (o12.score > o11.score and o12.score > min_score):
                    config_left = "flop"
                else:
                    config_left = "unclassified"
            
                if (o21.score > o22.score and o21.score > min_score):
                    config_right = "flip"
                elif (o22.score > o21.score and o22.score > min_score):
                    config_right = "flop"
                else:
                    config_right = "unclassified"
                    
            elif orientation == "right":
                forward_query = df[(df['flag'] == 0 ) | (df['flag'] == 2048)]['seq'].iloc[0]
                reverse_query = df[(df['flag'] == 16) | (df['flag'] == 2064)]['seq'].iloc[0]
                o11 = parasail.sw_trace(forward_query[-len(ff_seq.right_flip) - min_insert :], ff_seq.right_flip, 3, 1, SW_SCORE_MATRIX,)
                o12 = parasail.sw_trace(forward_query[-len(ff_seq.right_flop) - min_insert :], ff_seq.right_flop, 3, 1, SW_SCORE_MATRIX,)
                o21 = parasail.sw_trace(reverse_query[-len(ff_seq.right_flip) - min_insert :], ff_seq.right_flip, 3, 1, SW_SCORE_MATRIX,)
                o22 = parasail.sw_trace(reverse_query[-len(ff_seq.right_flop) - min_insert :], ff_seq.right_flop, 3, 1, SW_SCORE_MATRIX,)
    
                if o11.score > o12.score and o11.score > min_score:
                    config_right = "flip"
                elif (o12.score > o11.score and o12.score > min_score):
                    config_right = "flop"
                else:
                    config_right = "unclassified"
                if (o21.score > o22.score and o21.score > min_score):
                    config_left = "flip"
                elif (o22.score > o21.score and o22.score > min_score):
                    config_left = "flop"
                else:
                    config_left = "unclassified"  
                    
        elif len(df) == 1:
            query = df['seq'].iloc[0]
            if orientation == "left":    
                o11 = parasail.sw_trace(query, ff_seq.left_flip, 3, 1, SW_SCORE_MATRIX)
                o12 = parasail.sw_trace(query, ff_seq.left_flop, 3, 1, SW_SCORE_MATRIX)

                if o11.score > o12.score and o11.score > min_score:
                    config_left = "flip"
                elif (o12.score > o11.score and o12.score > min_score):
                    config_left = "flop"
                else:
                    config_left = "unclassified"
                
                config_right = "unclassified"
            
            elif orientation == "right":
                o11 = parasail.sw_trace(forward_query[-len(ff_seq.right_flip) - min_insert :], ff_seq.right_flip, 3, 1, SW_SCORE_MATRIX,)
                o12 = parasail.sw_trace(forward_query[-len(ff_seq.right_flop) - min_insert :], ff_seq.right_flop, 3, 1, SW_SCORE_MATRIX,)

                if o11.score > o12.score and o11.score > min_score:
                    config_right = "flip"
                elif (o12.score > o11.score and o12.score > min_score):
                    config_right = "flop"
                else:
                    config_right = "unclassified"

                config_left = "unclassified"
                
    elif (vector_type == "ss"):
        query = df['seq'].iloc[0]

        #if df["AX"].drop_duplicates().iloc[0] == "vector-partial":
        #if t["AX"] == "vector-partial":
            # ignore, since both sides are missing chunks of ITR
            #return "unclassified", "unclassified"

        if df["AX"].drop_duplicates().iloc[0] in ("vector-full", "vector-left-partial"):
            o1 = parasail.sw_trace(query, ff_seq.left_flip, 3, 1, SW_SCORE_MATRIX)
            o2 = parasail.sw_trace(query, ff_seq.left_flop, 3, 1, SW_SCORE_MATRIX)
            if o1.score > o2.score and o1.score > min_score:
                config_left = "flip"
            elif o2.score > o1.score and o2.score > min_score:
                config_left = "flop"
            else:
                config_left = "unclassified"

        if df["AX"].drop_duplicates().iloc[0] in ("vector-full", "vector-right-partial"):
            o1 = parasail.sw_trace(query[-len(ff_seq.right_flip) - 10 :],ff_seq.right_flip,3,1,SW_SCORE_MATRIX,)
            o2 = parasail.sw_trace(query[-len(ff_seq.right_flop) - 10 :],ff_seq.right_flop,3,1,SW_SCORE_MATRIX,)
            if o1.score > o2.score and o1.score > min_score:
                config_right = "flip"
            elif o2.score > o1.score and o2.score > min_score:
                config_right = "flop"
            else:
                config_right = "unclassified"
  
    return config_left, config_right


def load_per_read_info(fname):
    """Load per-read info, keyed by read IDs, from a CSV file."""
    with gzip.open(fname, "rt") as in_tsv:
        read_info = {r["read_id"]: r for r in csv.DictReader(in_tsv, delimiter="\t")}
    return read_info


def main(per_read_tsv, tagged_bam, vector_type, orientation, output_prefix, flipflop_fasta):
    """Entry point."""
    OUT_FIELDS = ["name", "type", "subtype", "start", "end", "leftITR", "rightITR"]

    if flipflop_fasta is None:
        flipflop_seqs = FlipFlopSeqSet(**SEQ_AAV2)
    else:
        flipflop_seqs = FlipFlopSeqSet.from_fasta(flipflop_fasta)

    read_info = load_per_read_info(per_read_tsv)

    with gzip.open(output_prefix + ".flipflop.tsv.gz", "wt") as fout:
            out_tsv = csv.writer(fout, delimiter="\t")
            out_tsv.writerow(OUT_FIELDS)
            reader = pysam.AlignmentFile(open(tagged_bam), "rb", check_sq=False)
            out_bam_full = pysam.AlignmentFile(
                open(output_prefix + ".flipflop-full.bam", "w"),
                "wb",
                header=reader.header,
            )
            out_bam_leftp = pysam.AlignmentFile(
                open(output_prefix + ".flipflop-left-partial.bam", "w"),
                "wb",
                header=reader.header,
            )
            out_bam_rightp = pysam.AlignmentFile(
                open(output_prefix + ".flipflop-right-partial.bam", "w"),
                "wb",
                header=reader.header,
            )

            # Initialize a list to hold the data
            data = []

            # Iterate through each record in reader
            for r in reader:
                d = r.to_dict()
                #print(r)
                d["header"] = r.header.to_dict()
                d["start"] = (r.reference_start)
                d["end"] = (r.reference_end)              
                data.append(d)

            # Convert the list of rows into a DataFrame
            df = pd.DataFrame(data)

            # Convert flag and tags columns to appropriate types
            df["flag"] = df["flag"].astype(int)
            df["tags"] = df["tags"].astype(str)
            
            #print(df)

            # Iterate through each read-id 'name' in the DataFrame
            for id in df['name'].drop_duplicates():
                print(id)
                f_df= df[df['name'] == id]

                # Iterate through each row of the DataFrame to extract AX and AT tags
                for index, row in f_df.iterrows():
                    tags = row["tags"]

                    # Extract AT value
                    at_match = re.search(r"AT:Z:(.+?)(?:,|'|$)", tags)
                    f_df.at[index, "AT"] = at_match.group(1) if at_match else None

                    # Extract AX value
                    ax_match = re.search(r"AX:Z:(.+?)(?:,|'|$)", tags)
                    f_df.at[index, "AX"] = ax_match.group(1) if ax_match else None

                if f_df["AT"].drop_duplicates().iloc[0] == "vector" and f_df["AX"].drop_duplicates().iloc[0] in (
                    "vector-full",
                    "vector-left-partial",
                    "vector-right-partial",
                ):
                    orientation= "left"
                    c_l, c_r = identify_flip_flop(f_df, flipflop_seqs, vector_type, orientation)
                    for index, row in f_df.iterrows():
                        if "AX" in f_df.columns and not f_df["AX"].empty:
                            if f_df["AX"].drop_duplicates().iloc[0] == "vector-full":
                                writer = out_bam_full
                            elif f_df["AX"].drop_duplicates().iloc[0] == "vector-right-partial":
                                writer = out_bam_leftp
                            elif f_df["AX"].drop_duplicates().iloc[0] == "vector-left-partial":
                                writer = out_bam_rightp
                        if row['name'] in read_info and "assigned_type" in read_info[row['name']]:
                            a_type = read_info[row['name']]["assigned_type"]
                        else:
                            a_type = None  # Handle missing or invalid data

                        f_df.at[index, "tags"] =   f_df.at[index, "tags"] + "AF:Z:" + c_l + "-" + c_r + "," + "AG:Z:" + a_type
                        if a_type not in ("scAAV", "ssAAV"):
                            continue

                        # Define the columns to filter
                        columns_to_keep = [
                            "name",
                            "flag",
                            "ref_name",
                            "ref_pos",
                            "map_quality",
                            "cigar",
                            "next_ref_name",
                            "next_ref_pos",
                            "length",
                            "seq",
                            "qual",
                            "tags",
                        ]

                        # Filter the DataFrame to keep only the specified columns
                        filtered_row = row[columns_to_keep]

                        header = pysam.AlignmentHeader.from_dict(row['header'])  # Convert OrderedDict to AlignmentHeader
                        filtered_row["flag"] = str(filtered_row["flag"])  # Convert directly if it's an integer

                        #r = pysam.AlignedSegment.from_dict(filtered_row, header)  # TO FIX             

                        out_tsv.writerow(
                            [
                                row['name'],  # Access qname from f_df
                                a_type,  # Access assigned_type from f_df
                                row["AX"][len("vector-") :],  # Extract AX value from f_df
                                str(row["start"]),  # Access reference_start from f_df
                                str(row["end"]),  # Access reference_end from f_df
                                c_l,  # Use c_l (computed earlier)
                                c_r,  # Use c_r (computed earlier)
                            ]
                        )
            

            out_bam_full.close()
            out_bam_leftp.close()
            out_bam_rightp.close()
            print("Output summmary:", fout.name)
    
    print(
        f"Individual BAM files written: {output_prefix}.vector- full,leftpartial,rightpartial -flipflop.bam"
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("sorted_tagged_bam", help="Sorted tagged BAM file")
    parser.add_argument("per_read_tsv", help="Per read TSV file")
    parser.add_argument("vector_type", help="AAV vector type")
    parser.add_argument("orientation", help="mITR orientation")
    parser.add_argument("-o", "--output-prefix", help="Output prefix", required=True)
    parser.add_argument(
        "--flipflop-fasta",
        default=None,
        help="(optional) flip flop fasta file (if not given, uses AAV2 default)",
    )

    args = parser.parse_args()

    main(
        args.per_read_tsv,
        args.sorted_tagged_bam,
        args.vector_type,
        args.orientation,
        args.output_prefix,
        args.flipflop_fasta,
    )
