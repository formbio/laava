#!/usr/bin/env python3
"""Get ITR flip flop configurations - Memory Optimized Version.

Must have already run `summarize_alignment.py` to get a .tagged.BAM file!
"""

from __future__ import annotations

import csv
import gzip
import gc
from typing import NamedTuple

import parasail
import pysam
from Bio import SeqIO
import re


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
        if not all(flipflip_seqs.values()):
            missing_keys = [
                key for key, val in known_seq_ids.items() if flipflip_seqs[val] is None
            ]
            raise RuntimeError(
                f"Required flip-flop sequence IDs missing from {fname}: " + ", ".join(
                    missing_keys
                )
            )
        return cls(**flipflip_seqs)


def ensure_bam_sorted_by_name(tagged_bam, output_prefix):
    """Ensure BAM file is sorted by read name."""
    try:
        # Check if BAM is already sorted by name
        header_lines = pysam.view("-H", tagged_bam).split("\n")
        is_sorted = any("SO:queryname" in line for line in header_lines)
        
        if not is_sorted:
            print("BAM file must be sorted by read name. Sorting now...")
            sorted_bam = output_prefix + ".name_sorted.bam"
            pysam.sort("-n", "-o", sorted_bam, tagged_bam)
            return sorted_bam
        
        return tagged_bam
    except Exception as e:
        print(f"Warning: Could not check BAM sort order ({e}). Proceeding with original file.")
        return tagged_bam


def is_vector_read(record):
    """Check if a record is a vector read that we need to process."""
    try:
        tags = dict(record.tags)
        at_tag = tags.get('AT', None)
        ax_tag = tags.get('AX', None)
        return (at_tag == 'vector' and 
                ax_tag in ('vector-full', 'vector-left-partial', 'vector-right-partial'))
    except:
        return False


def identify_flip_flop_direct(reads, ff_seq, vector_type, orientation):
    """Determine left and right flip/flop/unclassified configurations from pysam records directly."""
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
        if len(reads) == 2:
            if orientation == "left":    
                forward_read = next((r for r in reads if r.flag in (0, 2048)), None)
                reverse_read = next((r for r in reads if r.flag in (16, 2064)), None)
                
                if forward_read and reverse_read:
                    forward_query = forward_read.query_sequence
                    reverse_query = reverse_read.query_sequence
                    
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
                forward_read = next((r for r in reads if r.flag in (0, 2048)), None)
                reverse_read = next((r for r in reads if r.flag in (16, 2064)), None)
                
                if forward_read and reverse_read:
                    forward_query = forward_read.query_sequence
                    reverse_query = reverse_read.query_sequence
                    
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
                        
        elif len(reads) == 1:
            query = reads[0].query_sequence
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
                o11 = parasail.sw_trace(query[-len(ff_seq.right_flip) - min_insert :], ff_seq.right_flip, 3, 1, SW_SCORE_MATRIX,)
                o12 = parasail.sw_trace(query[-len(ff_seq.right_flop) - min_insert :], ff_seq.right_flop, 3, 1, SW_SCORE_MATRIX,)

                if o11.score > o12.score and o11.score > min_score:
                    config_right = "flip"
                elif (o12.score > o11.score and o12.score > min_score):
                    config_right = "flop"
                else:
                    config_right = "unclassified"

                config_left = "unclassified"
    
    # ssAAV support a single alignment read            
    elif (vector_type == "ss"):
        query = reads[0].query_sequence
        tags = dict(reads[0].tags)
        ax_tag = tags.get('AX', '')

        if ax_tag in ("vector-full", "vector-left-partial"):
            o1 = parasail.sw_trace(query, ff_seq.left_flip, 3, 1, SW_SCORE_MATRIX)
            o2 = parasail.sw_trace(query, ff_seq.left_flop, 3, 1, SW_SCORE_MATRIX)
            if o1.score > o2.score and o1.score > min_score:
                config_left = "flip"
            elif o2.score > o1.score and o2.score > min_score:
                config_left = "flop"
            else:
                config_left = "unclassified"

        if ax_tag in ("vector-full", "vector-right-partial"):
            o1 = parasail.sw_trace(query[-len(ff_seq.right_flip) - 10 :],ff_seq.right_flip,3,1,SW_SCORE_MATRIX,)
            o2 = parasail.sw_trace(query[-len(ff_seq.right_flop) - 10 :],ff_seq.right_flop,3,1,SW_SCORE_MATRIX,)
            if o1.score > o2.score and o1.score > min_score:
                config_right = "flip"
            elif o2.score > o1.score and o2.score > min_score:
                config_right = "flop"
            else:
                config_right = "unclassified"
  
    return config_left, config_right


def process_read_group(reads, ff_seq, vector_type, orientation, read_info, out_tsv, out_bam_full, out_bam_leftp, out_bam_rightp):
    """Process a group of 1-2 reads with the same name."""
    if not reads:
        return
    
    read_name = reads[0].query_name
    
    # Identify flip-flop configuration
    c_l, c_r = identify_flip_flop_direct(reads, ff_seq, vector_type, orientation)
    
    # Get assigned type and count
    a_type = None
    effective_count = 0
    if read_name in read_info and "assigned_type" in read_info[read_name]:
        a_type = read_info[read_name]["assigned_type"]
        effective_count = int(read_info[read_name]["effective_count"])
    
    # Process each read in the group
    for i, record in enumerate(reads):
        tags = dict(record.tags)
        ax_tag = tags.get('AX', '')
        
        # Add flip-flop tag
        new_tags = list(record.tags)
        new_tags.append(("AF", f"{c_l}-{c_r}"))
        if a_type:
            new_tags.append(("AG", a_type))
        record.tags = new_tags
        
        # Write to TSV if applicable (only once per group)
        if i == 0 and a_type in ("scAAV", "ssAAV"):
            for _ in range(effective_count):
                out_tsv.writerow([
                    read_name,
                    a_type,
                    ax_tag[len("vector-"):] if ax_tag.startswith("vector-") else ax_tag,
                    str(record.reference_start),
                    str(record.reference_end),
                    c_l,
                    c_r
                ])
        
        # Write to appropriate BAM file
        if ax_tag == "vector-full":
            out_bam_full.write(record)
        elif ax_tag == "vector-left-partial":
            out_bam_leftp.write(record)
        elif ax_tag == "vector-right-partial":
            out_bam_rightp.write(record)


# Function to load per-read information from a TSV file
def load_per_read_info(fname):
    """Load per-read info, keyed by read IDs, from a CSV file."""
    with gzip.open(fname, "rt") as in_tsv:
        read_info = {r["read_id"]: r for r in csv.DictReader(in_tsv, delimiter="\t")}
    return read_info


# Main function to process the tagged BAM file and per-read TSV file
def main(per_read_tsv, tagged_bam, vector_type, orientation, output_prefix, flipflop_fasta):
    """Entry point - Memory optimized streaming processing."""
    OUT_FIELDS = ["name", "type", "subtype", "start", "end", "leftITR", "rightITR"]

    if flipflop_fasta is None:
        flipflop_seqs = FlipFlopSeqSet(**SEQ_AAV2)
    else:
        flipflop_seqs = FlipFlopSeqSet.from_fasta(flipflop_fasta)

    # Ensure BAM is sorted by name
    sorted_bam = ensure_bam_sorted_by_name(tagged_bam, output_prefix)

    print("Loading per-read information...")
    read_info = load_per_read_info(per_read_tsv)
    print(f"Loaded {len(read_info):,} read records")

    with gzip.open(output_prefix + ".flipflop.tsv.gz", "wt") as fout:
        out_tsv = csv.writer(fout, delimiter="\t")
        out_tsv.writerow(OUT_FIELDS)
        
        with pysam.AlignmentFile(sorted_bam, "rb", check_sq=False) as reader:
            out_bam_full = pysam.AlignmentFile(
                output_prefix + ".flipflop-full.bam", "wb", header=reader.header)
            out_bam_leftp = pysam.AlignmentFile(
                output_prefix + ".flipflop-left-partial.bam", "wb", header=reader.header)
            out_bam_rightp = pysam.AlignmentFile(
                output_prefix + ".flipflop-right-partial.bam", "wb", header=reader.header)

            print("Starting memory-optimized streaming processing...")
            
            # Variables to track current processing state
            current_reads = []
            current_read_name = None
            total_processed = 0
            vector_reads_processed = 0
            
            for record in reader:
                total_processed += 1
                
                # Progress reporting
                if total_processed % 20000 == 0:
                    print(f"Processed {total_processed:,} total reads, {vector_reads_processed:,} vector reads")
                    break # for testing
                
                # Skip non-vector reads early to save memory
                if not is_vector_read(record):
                    continue
                
                read_name = record.query_name
                
                # Case 1: First read or continuing same group
                if current_read_name is None or read_name == current_read_name:
                    current_reads.append(record)
                    current_read_name = read_name
                    
                    # If we have 2 reads with same name, process them immediately
                    if len(current_reads) == 2:
                        process_read_group(current_reads, flipflop_seqs, vector_type, orientation, 
                                         read_info, out_tsv, out_bam_full, out_bam_leftp, out_bam_rightp)
                        vector_reads_processed += len(current_reads)
                        # Clear after processing
                        current_reads = []
                        current_read_name = None
                
                # Case 2: New read name encountered
                else:
                    # Process any pending reads from previous group
                    if current_reads:
                        process_read_group(current_reads, flipflop_seqs, vector_type, orientation, 
                                         read_info, out_tsv, out_bam_full, out_bam_leftp, out_bam_rightp)
                        vector_reads_processed += len(current_reads)
                    
                    # Start new group with this read
                    current_reads = [record]
                    current_read_name = read_name
            
            # Process any final pending reads
            if current_reads:
                process_read_group(current_reads, flipflop_seqs, vector_type, orientation, 
                                 read_info, out_tsv, out_bam_full, out_bam_leftp, out_bam_rightp)
                vector_reads_processed += len(current_reads)
            
            # Close BAM output files
            out_bam_full.close()
            out_bam_leftp.close()
            out_bam_rightp.close()
            
            print(f"Memory-optimized processing completed: {total_processed:,} total reads, {vector_reads_processed:,} vector reads processed")
            print("Output summary:", fout.name)
    
    print(
        f"Individual BAM files written: {output_prefix}.flipflop-full.bam, {output_prefix}.flipflop-left-partial.bam, {output_prefix}.flipflop-right-partial.bam"
    )
    
    # Clean up temporary sorted BAM if we created one
    if sorted_bam != tagged_bam:
        import os
        try:
            os.remove(sorted_bam)
            print(f"Cleaned up temporary sorted BAM: {sorted_bam}")
        except:
            print(f"Warning: Could not clean up temporary file: {sorted_bam}")


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
