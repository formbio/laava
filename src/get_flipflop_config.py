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
import pandas as pd


# Pre-compiled regex patterns for performance optimization
AT_PATTERN = re.compile(r"AT:Z:([^,'\s]+)")
AX_PATTERN = re.compile(r"AX:Z:([^,'\s]+)")

SW_SCORE_MATRIX = parasail.matrix_create("ACGT", 2, -5)

# Default AAV2 flip-flop sequences
SEQ_AAV2 = dict(
    left_flip="TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT",
    left_flop="TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT",
    right_flip="AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA",
    right_flop="AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA",
)

# Define a named tuple to hold the flip-flop sequences
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

def chunked_bam_reader(reader, chunk_size=10000):
    """Yield chunks of BAM records to prevent memory explosion on large datasets."""
    chunk = []
    for record in reader:
        chunk.append(record)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:  # Final partial chunk
        yield chunk


def process_dataframe_chunk(df, ff_seq, vector_type, orientation, read_info, out_tsv, out_bam_full, out_bam_leftp, out_bam_rightp):
    """Process a chunk of records using the existing DataFrame logic."""
    # Convert flag and tags columns to appropriate types
    df["flag"] = df["flag"].astype(int)
    df["tags"] = df["tags"].astype(str)
    
    # Extract tags using pre-compiled regex patterns for improved performance
    print(f"Extracting tags for {len(df)} records...")
    for index, row in df.iterrows():
        tags = row["tags"]
        # Use pre-compiled patterns for better performance than re.search()
        at_match = AT_PATTERN.search(tags)
        df.loc[index, "AT"] = at_match.group(1) if at_match else None
        ax_match = AX_PATTERN.search(tags)
        df.loc[index, "AX"] = ax_match.group(1) if ax_match else None
    
    # Filter to only vector reads to reduce processing overhead
    vector_reads = df[
        (df['AT'] == 'vector') & 
        (df['AX'].isin(['vector-full', 'vector-left-partial', 'vector-right-partial']))
    ]
    
    print(f"Processing {len(vector_reads)} vector reads from chunk...")
    
    # Use efficient groupby instead of quadratic operations
    for read_name, f_df in vector_reads.groupby('name'):
        if f_df["AT"].drop_duplicates().iloc[0] == "vector" and f_df["AX"].drop_duplicates().iloc[0] in (
            "vector-full",
            "vector-left-partial",
            "vector-right-partial",
        ):
            c_l, c_r = identify_flip_flop(f_df, ff_seq, vector_type, orientation)
            for index, row in f_df.iterrows():
                if row['name'] in read_info and "assigned_type" in read_info[row['name']]:
                    a_type = read_info[row['name']]["assigned_type"]
                    effective_count = int(read_info[row['name']]["effective_count"])
                else:
                    a_type = None  # Handle missing or invalid data
                    effective_count = 0

                f_df.at[index, "tags"] = f_df.at[index, "tags"] + "AF:Z:" + c_l + "-" + c_r + "," + "AG:Z:" + a_type
                
                if a_type not in ("scAAV", "ssAAV"):
                    continue

                for _ in range(effective_count):
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
                
                # Exit the loop after processing the first row
                break


# Function to identify flip-flop configurations based on the provided alingned reads and flip-flop sequences
def identify_flip_flop(df, ff_seq, vector_type, orientation):
    """Determine left and right flip/flip/unclassified configurations.

    Assume record tag:AT is vector, tag:AX can be full|left-partial|right-partial|partial
    Add back a tag 'AF' that is [flip/flop]-[flip/flop]
    """
    min_score = 250
    min_insert = 10 
    config_left, config_right = "unclassified", "unclassified"

    # Expected flags
    # 0 read forward strand
    # 16 read reverse strand
    # 2048 read forward strand + supplementary alignment
    # 2064 read reverse strand + supplementary alignment
    # 4 read unmapped
    
    # Full scAAV support two aligned reads, one for each side of the ITR one for truncation before the mITR
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
    
    # ssAAV support a single alignment read            
    elif (vector_type == "ss"):
        query = df['seq'].iloc[0]

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

# Function to load per-read information from a TSV file
def load_per_read_info(fname):
    """Load per-read info, keyed by read IDs, from a CSV file."""
    with gzip.open(fname, "rt") as in_tsv:
        read_info = {r["read_id"]: r for r in csv.DictReader(in_tsv, delimiter="\t")}
    return read_info

# Main function to process the tagged BAM file and per-read TSV file
def main(per_read_tsv, tagged_bam, vector_type, orientation, output_prefix, flipflop_fasta):
    """Entry point - Chunked processing to prevent OOM on large datasets."""
    OUT_FIELDS = ["name", "type", "subtype", "start", "end", "leftITR", "rightITR"]
    CHUNK_SIZE = 20000  # Process 10K records at a time to prevent memory explosion

    if flipflop_fasta is None:
        flipflop_seqs = FlipFlopSeqSet(**SEQ_AAV2)
    else:
        flipflop_seqs = FlipFlopSeqSet.from_fasta(flipflop_fasta)

    print("Loading per-read information...")
    read_info = load_per_read_info(per_read_tsv)
    print(f"Loaded {len(read_info):,} read records")

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

        print(f"Starting chunked processing (chunk size: {CHUNK_SIZE:,} records)...")
        chunk_count = 0
        total_processed = 0
        
        # CHUNKED PROCESSING: Process BAM in manageable chunks to prevent OOM
        for chunk_records in chunked_bam_reader(reader, CHUNK_SIZE):
            chunk_count += 1
            chunk_size = len(chunk_records)
            total_processed += chunk_size
            
            print(f"Processing chunk {chunk_count} ({chunk_size:,} records, {total_processed:,} total)...")
            
            # Convert chunk to DataFrame (much smaller than full BAM)
            chunk_data = []
            for r in chunk_records:
                d = r.to_dict()
                d["header"] = r.header.to_dict()
                d["start"] = r.reference_start
                d["end"] = r.reference_end              
                chunk_data.append(d)
            
            # Create DataFrame from chunk (limited memory usage)
            df_chunk = pd.DataFrame(chunk_data)
            
            # Process chunk using existing logic
            process_dataframe_chunk(
                df_chunk, flipflop_seqs, vector_type, orientation, 
                read_info, out_tsv, out_bam_full, out_bam_leftp, out_bam_rightp
            )
            break # only do one chunk for testing

        out_bam_full.close()
        out_bam_leftp.close()
        out_bam_rightp.close()
        print(f"Chunked processing completed: {chunk_count} chunks, {total_processed:,} total records")
        print("Output summary:", fout.name)
    
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
