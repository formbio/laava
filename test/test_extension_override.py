#!/usr/bin/env python3
"""
Simple Extension Override for Testing

This demonstrates the configuration-driven function pointer system with a minimal,
safe override that adds logging to show the extension system is working.

This is identical to the original iter_cigar_w_aligned_pair function but with
clear logging messages to prove the extension was loaded and executed.
"""

import logging
import itertools

# CIGAR operation mapping (copied from summarize_alignment.py)
CIGAR_DICT = {
    0: "M",
    1: "I", 
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    7: "=",
    8: "X",
}

def iter_cigar(rec):
    """Iterator for CIGAR operations (copied from summarize_alignment.py)."""
    for intcode, count in rec.cigartuples:
        if intcode not in CIGAR_DICT:
            raise RuntimeError(f"Unexpected cigar string: {rec.cigarstring}")
        charcode = CIGAR_DICT[intcode]
        if charcode == "H":
            # Ignore hard-clipping (due to supp alignment)
            continue
        yield from itertools.repeat((charcode, count), count)

def iter_cigar_w_aligned_pair_with_logging(rec, writer):
    """
    Extension override with logging to demonstrate function pointer system.
    
    This function is IDENTICAL to the original iter_cigar_w_aligned_pair
    but adds clear logging messages to prove the extension system works.
    
    Args:
        rec: pysam AlignedSegment record
        writer: CSV writer for mismatch output
        
    Returns:
        Tuple of (total_errors, total_length)
    """
    # CLEAR LOGGING TO SHOW EXTENSION IS ACTIVE
    logging.info("🔧 START of override extension! (test_extension_override.py)")
    
    # IDENTICAL LOGIC TO ORIGINAL FUNCTION
    prev_cigar_type = None
    prev_r_pos = 0
    total_err = 0
    total_len = 0
    
    # ENH: itertools.zip_longest + safety check
    for (_q_pos, r_pos), (cigar_type, cigar_count) in zip(
        rec.get_aligned_pairs(), iter_cigar(rec)
    ):
        if cigar_type == "S":
            # Nothing to do if soft-clipped, r_pos must be None
            if r_pos is not None:
                logging.warning("Unexpected value for r_pos: %s", r_pos)
            continue
        total_len += cigar_count
        if cigar_type != prev_cigar_type:
            if cigar_type in ("I", "D", "X", "N"):
                total_err += cigar_count
                info = {
                    "read_id": rec.qname,
                    "pos0": r_pos if cigar_type != "I" else prev_r_pos,
                    "type": cigar_type,
                    "type_len": cigar_count,
                }
                writer.writerow(info)
            prev_cigar_type = cigar_type
        if r_pos is not None:
            prev_r_pos = r_pos
    
    # CLEAR LOGGING TO SHOW EXTENSION COMPLETED
    logging.info("🔧 END of override extension! (test_extension_override.py)")
    
    return total_err, total_len

def get_extension_info():
    """
    Get information about this test extension.
    
    Returns:
        Dictionary with extension details
    """
    return {
        "name": "Test Extension Override",
        "version": "1.0.0", 
        "description": "Simple logging override to demonstrate extension system",
        "purpose": "Testing and validation",
        "identical_to_original": True,
        "adds_logging": True
    }
