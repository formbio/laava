#!/usr/bin/env python3
"""Summarize the AAV alignment."""

import gzip
import itertools
import logging
import os
import re
import shutil
import subprocess
import sys
from csv import DictReader, DictWriter
from multiprocessing import Process
# import pdb

import pysam

CIGAR_DICT = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X",
    9: "B",
}

annot_rex = re.compile(r"NAME=(\S+);TYPE=([a-zA-Z]+);(REGION=\d+\-\d+){0,1}")
ccs_rex = re.compile(r"\S+\/\d+\/ccs(\/fwd|\/rev)?")
ANNOT_TYPE_PRIORITIES = {"vector": 1, "repcap": 2, "helper": 3, "lambda": 4, "host": 5}

MAX_MISSING_FLANK = 100
MAX_OUTSIDE_VECTOR = 100
TARGET_GAP_THRESHOLD = 200  # skipping through the on-target region for more than this is considered "full-gap"


def subset_sam_by_readname_list(
    in_bam,
    out_bam,
    per_read_tsv,
    wanted_types,
    wanted_subtypes,
    max_count=None,
    exclude_subtype=False,
    exclude_type=False,
):
    qname_list = {}  # qname --> (a_type, a_subtype)
    for r in DictReader(open(per_read_tsv), delimiter="\t"):
        # pdb.set_trace()
        if (
            wanted_types is None
            or (not exclude_type and r["assigned_type"] in wanted_types)
            or (exclude_type and r["assigned_type"] not in wanted_types)
        ) and (
            wanted_subtypes is None
            or (not exclude_subtype and (r["assigned_subtype"] in wanted_subtypes))
            or (exclude_subtype and (r["assigned_subtype"] not in wanted_subtypes))
        ):
            qname_list[r["read_id"]] = (r["assigned_type"], r["assigned_subtype"])

    cur_count = 0
    reader = pysam.AlignmentFile(in_bam, "rb", check_sq=False)
    writer = pysam.AlignmentFile(out_bam, "wb", header=reader.header)
    for r in reader:
        if r.qname in qname_list:
            d = add_assigned_types_to_record(r, *qname_list[r.qname])
            writer.write(pysam.AlignedSegment.from_dict(d, reader.header))
            cur_count += 1
            if max_count is not None and cur_count >= max_count:
                break
    reader.close()
    writer.close()


def iter_cigar(rec):
    # first we exclude cigar front/end that is hard-clipped (due to supp alignment)
    cigar_list = rec.cigar
    if CIGAR_DICT[cigar_list[0][0]] == "H":
        cigar_list = cigar_list[1:]
    if CIGAR_DICT[cigar_list[-1][0]] == "H":
        cigar_list = cigar_list[:-1]

    # warning_printed = False
    for _type, count in cigar_list:
        x = CIGAR_DICT[_type]
        if x in ("M", "=", "X", "I", "D", "S", "N"):
            yield from itertools.repeat((x, count), count)
        else:
            raise RuntimeError(f"Unexpected cigar {count}{x} seen! Abort!")


def iter_cigar_w_aligned_pair(rec, writer):
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
    return total_err, total_len


def read_annotation_file(annot_filename):
    """Read the annotation.txt file into a dictionary.

    Example:

    NAME=chr1;TYPE=host;
    NAME=chr2;TYPE=host;
    NAME=myVector;TYPE=vector;REGION=1795-6553;
    NAME=myCapRep;TYPE=repcap;REGION=1895-5987;
    NAME=myHelper;TYPE=helper;

    :param annot_filename: Annotation file following the format indicated above. Only "vector" is required. Others optional.
    :return:
    """
    result = {}
    for line in open(annot_filename):
        stuff = line.strip()
        m = annot_rex.match(stuff)
        if m is None:
            raise RuntimeError(
                f"{stuff} is not a valid annotation line! Should follow format "
                "`NAME=xxxx;TYPE=xxxx;REGION=xxxx;`. Abort!"
            )

        seq_name = m.group(1)
        ref_type = m.group(2)
        coord_region = (
            None
            if m.group(3) is None
            else tuple(map(int, m.group(3).split("=")[1].split("-")))
        )
        if ref_type in result:
            raise RuntimeError(f"Annotation file has multiple {ref_type} types. Abort!")
        if ref_type not in ANNOT_TYPE_PRIORITIES:
            raise RuntimeError(
                f"{ref_type} is not a valid type (host, repcap, vector, helper). Abort!"
            )
        result[seq_name] = {"type": ref_type, "region": coord_region}
    return result


def is_on_target(r, target_start, target_end):
    """Determine the alignment subtype.

    Possible assign types:
     - full (within target_start, target_end)
     - backbone (outside target_start, target_end)
     - left-partial (incomplete only on the 3'/right end)
     - right-partial (incomplete only on the 5'/left end)
     - partial (incomplete both on the 5' and 3' end)
     - vector+backbone (any amount that crosses over on target and backbone)
    """
    # Too far beyond left/5' target start
    left_long = r.reference_start < (target_start - MAX_OUTSIDE_VECTOR)
    # Incomplete on left/5'
    left_short = (target_start + MAX_MISSING_FLANK) < r.reference_start
    # Complete 5' start/left
    left_ok = not left_short and not left_long

    # Incomplete on right/3'
    right_short = r.reference_end < (target_end - MAX_MISSING_FLANK)
    # Too far beyond right/3' target end
    right_long = (target_end + MAX_OUTSIDE_VECTOR) < r.reference_end
    # Complete 3' end/right
    right_ok = not right_short and not right_long

    # Check the common case -- in most samples, most reads are on target
    if left_ok:
        if right_ok:
            for cigar_type, num in iter_cigar(r):
                if cigar_type == "N" and num >= TARGET_GAP_THRESHOLD:
                    return "full-gap"
            return "full"
        if right_short:
            return "left-partial"
        # Into backbone on right/3'
        if not right_long:
            logging.warning("Unexpected value for right_long: %s", right_long)
        return "vector+backbone"
    # Check for no overlap with target
    if r.reference_end < target_start or r.reference_start > target_end:
        return "backbone"
    # Check remaining cases
    if left_short:
        # Incomplete left/5' -> look at right for partial or right-partial
        if right_ok:
            return "right-partial"
        if right_short:
            return "partial"
    # Overlap both target and backbone, at either end
    if not (
        (r.reference_start < target_start < r.reference_end)
        or (r.reference_start < target_end < r.reference_end)
    ):
        logging.warning(
            "Unexpected vector+backbone alignment vs. target coordinates: %s vs. %s",
            (r.reference_start, r.reference_end),
            (target_start, target_end),
        )
    return "vector+backbone"


def assign_read_type(r, annotation):
    """Determine the read alignment type and subtype.

    :param read_dict: dict of {'supp', 'primary'}
    :return: assigned_type, which could be (scAAV, ssAAV, unclassified) + (super, full, partial, unclassified)

    <assigned_type: ssAAV, scAAV, backbone, helper, repcap, host, can use “+” sign>,
    <assigned_subtype: full, partial, nonAAV>
    <map_stat: unmapped | fully_aligned | partial_aligned | chimeric_aligned>,
    <map to: comma-delimited list of [chr:start-end]>,
    <comma-delimited list of unmapped portion, if any>,
    """
    read_type = annotation[r.reference_name]["type"]
    if annotation[r.reference_name]["region"] is not None:
        return read_type, is_on_target(r, *annotation[r.reference_name]["region"])
    else:
        return read_type, "NA"


def process_alignment_bam(
    sample_id,
    sorted_sam_filename,
    annotation,
    output_prefix,
    starting_readname=None,
    ending_readname=None,
):
    """Process the read alignments versus annotations.

    :param sorted_sam_filename: Sorted (by read name) SAM filename
    :param annotation:
    :param output_prefix:
    """
    ALIGNMENT_FIELDS = [
        "sample_id",
        "read_id",
        "read_len",
        "is_mapped",
        "is_supp",
        "map_name",
        "map_start0",
        "map_end1",
        "map_len",
        "map_iden",
        "map_type",
        "map_subtype",
    ]
    PER_READ_FIELDS = [
        "read_id",
        "read_len",
        "has_primary",
        "has_supp",
        "assigned_type",
        "assigned_subtype",
        "effective_count",
    ]
    NONMATCH_FIELDS = ["read_id", "pos0", "type", "type_len"]

    f_alignments = open(output_prefix + ".alignments.tsv", "w")
    f_nonmatch = open(output_prefix + ".nonmatch_stat.tsv", "w")
    f_per_read = open(output_prefix + ".per_read.tsv", "w")

    out_alignments = DictWriter(f_alignments, ALIGNMENT_FIELDS, delimiter="\t")
    out_nonmatch = DictWriter(f_nonmatch, NONMATCH_FIELDS, delimiter="\t")
    out_per_read = DictWriter(f_per_read, PER_READ_FIELDS, delimiter="\t")
    out_alignments.writeheader()
    out_nonmatch.writeheader()
    out_per_read.writeheader()

    reader = pysam.AlignmentFile(sorted_sam_filename, check_sq=False)
    bam_writer = pysam.AlignmentFile(
        output_prefix + ".tagged.bam", "wb", header=reader.header
    )

    if starting_readname is not None:
        # progress forward until we get to the read
        while True:
            try:
                cur_r = next(reader)
                if cur_r.qname == starting_readname:
                    records = [cur_r]
                    break
            except StopIteration:
                break
    else:
        records = [
            next(reader)
        ]  # records will hold all the multiple alignment records of the same read

    while True:
        try:
            cur_r = next(reader)
            if ending_readname is not None and cur_r.qname == ending_readname:
                break
            if cur_r.qname != records[-1].qname:
                process_alignment_records_for_a_read(
                    sample_id,
                    records,
                    annotation,
                    out_alignments,
                    out_nonmatch,
                    out_per_read,
                    bam_writer,
                )
                records = [cur_r]
            else:
                records.append(cur_r)
        except StopIteration:  # finished reading the SAM file
            break

    # Finish processing the last records
    process_alignment_records_for_a_read(
        sample_id,
        records,
        annotation,
        out_alignments,
        out_nonmatch,
        out_per_read,
        bam_writer,
    )
    bam_writer.close()
    f_alignments.close()
    f_nonmatch.close()
    f_per_read.close()
    return f_per_read.name, output_prefix + ".tagged.bam"


MIN_PRIM_SUPP_COV = 0.8  # at minimum the total of prim + main supp should cover this much of the original sequence


def find_companion_supp_to_primary(prim, supps):
    """Return the most likely companion supp to the primary.

    :param prim: the primary info
    :param supps: the list of supp info
    :return: return the most likely companion supp to the primary
    """

    def get_true_start_end(rec, true_qlen):
        # first we need to look at the strand
        # then also look at clipping
        cigartype, cigarlen = rec.cigartuples[0]
        offset = cigarlen if CIGAR_DICT[cigartype] == "H" else 0
        if rec.is_reverse:  # on - strand
            # we need to know the true length
            return true_qlen - (rec.qend + offset), true_qlen - (rec.qstart + offset)
        else:  # on + strand
            # just need to look at clipping
            return rec.qstart + offset, rec.qend + offset

    # if prim['rec'].qname=='m64011_220616_211638/9503552/ccsfwd':
    # pdb.set_trace()
    supp_should_be_rev = not prim["rec"].is_reverse
    # first look for a +/- supp
    for supp in supps:
        if supp["rec"].is_reverse == supp_should_be_rev:
            prim_start, prim_end = get_true_start_end(prim["rec"], prim["read_len"])
            supp_start, supp_end = get_true_start_end(supp["rec"], supp["read_len"])
            min_start = min(prim_start, supp_start)
            max_end = max(prim_end, supp_end)
            if (max_end - min_start) >= prim["read_len"] * MIN_PRIM_SUPP_COV:
                return supp, "+/-"

    # if that didn't work, check if there's a +/+ supp
    for supp in supps:
        if supp["rec"].is_reverse == prim["rec"].is_reverse:
            prim_start, prim_end = get_true_start_end(prim["rec"], prim["read_len"])
            supp_start, supp_end = get_true_start_end(supp["rec"], supp["read_len"])
            min_start = min(prim_start, supp_start)
            max_end = max(prim_end, supp_end)
            if (max_end - min_start) >= prim["read_len"] * MIN_PRIM_SUPP_COV:
                return supp, "+/+"
    return None, None


def add_assigned_types_to_record(r, a_type, a_subtype):
    """Add BAM tags.

    AT tag <type:scAAV|ssAAV|unclassified>
    AS tag <type:>
    AX tag which is "AT-AX"
    """
    d = r.to_dict()
    d["tags"].append("AT:Z:" + a_type)
    d["tags"].append("AS:Z:" + a_subtype)
    d["tags"].append("AX:Z:" + a_type + "-" + a_subtype)
    return d


def process_alignment_records_for_a_read(
    sample_id,
    records,
    annotation,
    out_alignments,
    out_nonmatch,
    out_per_read,
    bam_writer,
):
    """For each, find the most probable assignment.

    Priority: vector > rep/cap > helper > host

    :param records: list of alignment records for the same read
    :return:
    """
    read_tally = {"primary": None, "supp": []}
    for r in records:
        # check ccs id format is <movie>/<zmw>/ccs[/rev or /fwd]
        if ccs_rex.fullmatch(r.qname) is None:
            logging.warning(
                "WARNING: sequence ID does not follow format movie/zmw/ccs[/rev or /fwd]. "
                "Might undercount ssAAV!"
            )

        info = {
            "sample_id": sample_id,
            "read_id": r.qname,
            "read_len": r.query_length,
            "is_mapped": "N" if r.is_unmapped else "Y",
            "is_supp": "NA",
            "rec": r,  # we won't write this out later, it's a holder here for processing prim v supp
            "map_name": "NA",
            "map_start0": "NA",
            "map_end1": "NA",
            "map_len": "NA",
            "map_iden": "NA",
            "map_type": "NA",
            "map_subtype": "NA",
        }
        if r.is_unmapped:
            read_tally["primary"] = info
        else:
            cigar_list = r.cigar
            seq_len = r.query_length
            if CIGAR_DICT[cigar_list[0][0]] == "H":
                seq_len += cigar_list[0][1]
            if CIGAR_DICT[cigar_list[-1][0]] == "H":
                seq_len += cigar_list[-1][1]

            info["map_name"] = r.reference_name
            info["read_len"] = seq_len
            info["is_supp"] = "Y" if r.is_supplementary else "N"
            info["map_start0"] = r.reference_start
            info["map_end1"] = r.reference_end
            info["map_len"] = r.reference_end - r.reference_start
            total_err, total_len = iter_cigar_w_aligned_pair(r, out_nonmatch)
            info["map_iden"] = format(1.0 - total_err / total_len, ".10f")

            a_type, a_subtype = assign_read_type(r, annotation)
            info["map_type"] = a_type
            info["map_subtype"] = a_subtype
            logging.debug("%s %s %s", r.qname, a_type, a_subtype)
            # pdb.set_trace()

            if r.is_supplementary:
                read_tally["supp"].append(info)
            else:
                assert read_tally["primary"] is None
                read_tally["primary"] = info
        # NB: will write primary/supp to `out_alignments` after ruling out non-compatible subs

    # summarize it per read, now that all relevant alignments have been processed
    prim = read_tally["primary"]
    supps = read_tally["supp"]

    if len(supps) == 0:
        supp = None
        supp_orientation = None
    elif len(supps) >= 1:  # there's multiple supp, find the companion matching supp
        supp, supp_orientation = find_companion_supp_to_primary(prim, supps)
        # supp could be None, in which case there is best matching supp!
        # in the case supp is None we wanna see if this is a weird read (ex: mapped twice to + strand)

    # write the assigned type / subtype to the new BAM output
    bam_writer.write(
        pysam.AlignedSegment.from_dict(
            add_assigned_types_to_record(
                prim["rec"], prim["map_type"], prim["map_subtype"]
            ),
            prim["rec"].header,
        )
    )
    del prim["rec"]
    out_alignments.writerow(prim)
    if supp is not None:
        bam_writer.write(
            pysam.AlignedSegment.from_dict(
                add_assigned_types_to_record(
                    supp["rec"], supp["map_type"], supp["map_subtype"]
                ),
                supp["rec"].header,
            )
        )
        del supp["rec"]
        out_alignments.writerow(supp)

    sum_info = {
        "read_id": prim["read_id"],
        "read_len": prim["read_len"],
        "has_primary": prim["is_mapped"],
        "has_supp": "Y" if supp is not None else "N",
        "assigned_type": "NA",
        "assigned_subtype": "NA",
        "effective_count": 1,
    }
    if sum_info["has_primary"] == "Y":
        if prim["map_type"] == "vector":
            if supp is None:  # special case: primary only, maps to vector --> is ssAAV
                # double check the special case where there was supp candidates but no companion
                if len(supps) > 0:
                    sum_type = "unclassified"  # might be a weird case ex: a read covers the region twice as on + strand
                    sum_subtype = prim["map_subtype"]
                else:  # never had any supp candidates, def ssAAV
                    sum_type = "ssAAV"
                    sum_subtype = prim["map_subtype"]
            elif supp["map_type"] == "vector":
                if supp_orientation == "+/-":
                    # special case, primary+ supp, maps to vector, --> is scAAV
                    sum_type = "scAAV"
                else:
                    assert supp_orientation == "+/+"
                    sum_type = "tandem"
                if supp["map_subtype"] == prim["map_subtype"]:
                    sum_subtype = prim["map_subtype"]
                else:
                    sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
            else:  # primary is in vector, supp not in vector
                sum_type = prim["map_type"] + "|" + supp["map_type"]
                sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
        # mapping to non-AAV vector region
        elif supp is None:
            sum_type = prim["map_type"]
            sum_subtype = prim["map_subtype"]
        elif supp["map_type"] == prim["map_type"]:
            sum_type = prim["map_type"]
            if supp["map_subtype"] == prim["map_subtype"]:
                sum_subtype = prim["map_subtype"]
            else:
                sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
        else:
            sum_type = prim["map_type"] + "|" + supp["map_type"]
            sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
        sum_info["assigned_type"] = sum_type
        sum_info["assigned_subtype"] = sum_subtype

    # ###############################################################
    # 'effective_count' - now look at whether this is an ssAAV type
    # ex: <movie>/<zmw>/ccs   means potential two species (effective count of 2)
    # ex: <movie>/<zmw>/ccs/fwd or rev   is just one
    # NOTE: this can still be undercounting becuz not considering thing that are not ssAAV
    # ###############################################################
    if sum_info["assigned_type"] == "ssAAV":
        if sum_info["read_id"].endswith("/ccs"):
            sum_info["effective_count"] = 2
        # elif sum_info['read_id'].endswith('/ccs/fwd') or sum_info['read_id'].endswith('/ccs/rev'):
        #    sum_info['effective_count'] = 1  # not needed, default is to 1

    out_per_read.writerow(sum_info)
    logging.debug("%s", sum_info)
    # pdb.set_trace()


def run_processing_parallel(
    sample_id, sorted_sam_filename, annotation, output_prefix, num_chunks=1
):
    reader = pysam.AlignmentFile(open(sorted_sam_filename), check_sq=False)
    readname_list = [next(reader).qname]
    for r in reader:
        if r.qname != readname_list[-1]:
            readname_list.append(r.qname)

    total_num_reads = len(readname_list)
    chunk_size = (total_num_reads // num_chunks) + 1
    logging.info(
        f"Total {total_num_reads} reads, dividing into {num_chunks} "
        f"chunks of size {chunk_size}..."
    )

    pool = []
    for i in range(num_chunks):
        p = Process(
            target=process_alignment_bam,
            args=(
                sample_id,
                sorted_sam_filename,
                annotation,
                output_prefix + "." + str(i + 1),
                readname_list[i * chunk_size],
                None
                if (i + 1) * chunk_size > total_num_reads
                else readname_list[(i + 1) * chunk_size],
            ),
        )
        p.start()
        pool.append(p)
        logging.info("Going from %s to %s", i * chunk_size, (i + 1) * chunk_size)
    for i, p in enumerate(pool):
        logging.debug(f"DEBUG: Waiting for {i}th pool to finish.")
        p.join()

    # Combine the data together for *.nonmatch_stat.tsv, *.per_read.tsv, *.alignments.tsv

    # Closes over: output_prefix, num_chunks
    def gather_text_chunks(suffix, compress=False):
        logging.info("Combining chunk data... (*%s)", suffix)
        # Copy the first chunk over
        if compress:
            out_path = output_prefix + suffix + ".gz"
            f_out = gzip.open(out_path, "wb")
        else:
            out_path = output_prefix + suffix
            f_out = open(out_path, "w")
        first_chunk = f"{output_prefix}.1{suffix}"
        chunk_paths = [first_chunk]
        with open(first_chunk) as f_in:
            if compress:
                for line in f_in:
                    f_out.write(line.encode())
            else:
                shutil.copyfileobj(f_in, f_out)
        # Copy the remaining chunks
        for i in range(1, num_chunks):
            chunk_path = f"{output_prefix}.{i+1}{suffix}"
            with open(chunk_path) as f_in:
                f_in.readline()  # Skip the header
                if compress:
                    for line in f_in:
                        f_out.write(line.encode())
                else:
                    shutil.copyfileobj(f_in, f_out)
            chunk_paths.append(chunk_path)
        f_out.close()
        # Delete the chunk data
        logging.info("Data combining complete. Deleting chunk data (*%s).", suffix)
        for chunk_path in chunk_paths:
            os.remove(chunk_path)
        return out_path

    _outpath_nonmatch = gather_text_chunks(
        ".nonmatch_stat.tsv", compress=True
    )  # -> .nonmatch_stat.tsv.gz
    outpath_per_read = gather_text_chunks(".per_read.tsv")
    _outpath_alignments = gather_text_chunks(".alignments.tsv")

    # Combine the data together for *.tagged.bam
    # TODO - samtools cat/merge?
    logging.info("Combining BAM chunk data...")
    # Copy the first chunk over
    first_bam_chunk = f"{output_prefix}.1.tagged.bam"
    bam_chunk_paths = [first_bam_chunk]
    bam_reader = pysam.AlignmentFile(first_bam_chunk, "rb", check_sq=False)
    outpath_bam = output_prefix + ".tagged.bam"
    f_tagged_bam = pysam.AlignmentFile(outpath_bam, "wb", template=bam_reader)
    for r in bam_reader:
        f_tagged_bam.write(r)
    # Copy the remaining chunks
    for i in range(1, num_chunks):
        chunk_path = f"{output_prefix}.{i+1}.tagged.bam"
        for r in pysam.AlignmentFile(chunk_path, "rb", check_sq=False):
            f_tagged_bam.write(r)
        bam_chunk_paths.append(chunk_path)
    f_tagged_bam.close()
    bam_reader.close()
    # Delete the chunk data
    logging.info("Data combining complete. Deleting BAM chunk data.")
    for chunk_path in bam_chunk_paths:
        os.remove(chunk_path)

    return outpath_per_read, outpath_bam


def main(args):
    """Entry point."""
    annotation = read_annotation_file(args.annotation_txt)
    if args.cpus == 1:
        per_read_tsv, full_out_bam = process_alignment_bam(
            args.sample_id, args.sam_filename, annotation, args.output_prefix
        )
    else:
        per_read_tsv, full_out_bam = run_processing_parallel(
            args.sample_id,
            args.sam_filename,
            annotation,
            args.output_prefix,
            num_chunks=args.cpus,
        )

    # subset BAM files into major categories for ease of loading into IGV for viewing
    # subset_sam_by_readname_list(in_bam, out_bam, per_read_tsv, wanted_types, wanted_subtypes)
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".scAAV-full.tagged.bam",
        per_read_tsv,
        ["scAAV"],
        ["full"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".scAAV-partials.tagged.bam",
        per_read_tsv,
        ["scAAV"],
        ["partial", "left-partial", "right-partial"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".scAAV-other.tagged.bam",
        per_read_tsv,
        ["scAAV"],
        ["partial", "left-partial", "right-partial", "full"],
        exclude_subtype=True,
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".ssAAV-full.tagged.bam",
        per_read_tsv,
        ["ssAAV"],
        ["full"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".ssAAV-partials.tagged.bam",
        per_read_tsv,
        ["ssAAV"],
        ["partial", "left-partial", "right-partial"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".ssAAV-other.tagged.bam",
        per_read_tsv,
        ["ssAAV"],
        ["partial", "left-partial", "right-partial", "full"],
        exclude_subtype=True,
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".others.tagged.bam",
        per_read_tsv,
        ["ssAAV", "scAAV"],
        None,
        exclude_type=True,
    )

    # samtools sort/index the above files
    try:
        subprocess.check_call("samtools --help > /dev/null", shell=True)
    except Exception:
        logging.warning(
            "WARNING: unable to call samtools to sort the output BAM files. End."
        )
        sys.exit(-1)

    parts = [
        ".scAAV-full",
        ".scAAV-partials",
        ".scAAV-other",
        ".ssAAV-full",
        ".ssAAV-partials",
        ".ssAAV-other",
        ".others",
    ]
    for part in parts:
        p = args.output_prefix + part
        subprocess.check_call(
            f"samtools sort {p}.tagged.bam > {p}.tagged.sorted.bam", shell=True
        )
        subprocess.check_call(f"samtools index {p}.tagged.sorted.bam", shell=True)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("sam_filename", help="Sorted by read name SAM file")
    parser.add_argument("annotation_txt", help="Annotation file")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("-i", "--sample-id", required=True, help="Sample unique ID")
    parser.add_argument(
        "--max-allowed-missing-flanking",
        default=100,
        type=int,
        help="""Maximum allowed, in bp, missing from the left and right flanks of the
                annotated vector region to be still considered 'full' rather than
                'partial'.
                [Default: %(default)s]""",
    )
    parser.add_argument(
        "--max-allowed-outside-vector",
        default=100,
        type=int,
        help="""Maximum allowed, in bp, outside the annotated vector region for a read
                that fully covers the vector region to be still considered
                'full' rather than 'vector+backbone'.
                [Default: %(default)s]""",
    )
    parser.add_argument(
        "--target-gap-threshold",
        default=200,
        type=int,
        help="""Skipping through the on-target region for more than this is
                considered 'full-gap'. [Default: %(default)s]""",
    )
    parser.add_argument(
        "--cpus", default=1, type=int, help="Number of CPUs. [Default: %(default)s]"
    )
    parser.add_argument("--debug", action="store_true", default=False)

    args = parser.parse_args()

    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=log_level, format="%(message)s")
    MAX_MISSING_FLANK = args.max_allowed_missing_flanking
    MAX_OUTSIDE_VECTOR = args.max_allowed_outside_vector
    TARGET_GAP_THRESHOLD = args.target_gap_threshold

    main(args)
