#!/usr/bin/env python3
"""Print target region start/end coords for vector and repcap regions."""

from prepare_annotation import read_annotation_bed


if __name__ == "__main__":
    from argparse import ArgumentParser

    AP = ArgumentParser(description=__doc__)
    AP.add_argument("annotation_bed", help="Annotation file")
    AP.add_argument("itr_labels", nargs="*", help="ITR label(s) in annotation BED")
    args = AP.parse_args()

    ann = read_annotation_bed(args.annotation_bed, args.itr_labels)
    vec_start1 = ann["vector"].start1
    vec_end = ann["vector"].end
    if ann["repcap"] is None:
        rc_start1 = 0
        rc_end = 0
    else:
        rc_start1 = ann["repcap"].start1
        rc_end = ann["repcap"].end
    print(vec_start1, vec_end, rc_start1, rc_end)
