#!/usr/bin/env python3
"""Convert legacy annotation.txt to UCSC BED."""

import argparse
import csv
import logging
import re

RX_ANNOT = re.compile(r"NAME=([^;\s]+);TYPE=([^;\s]+);(REGION=\d+\-\d+){0,1}")
ANNOT_TYPE_PRIORITIES = {"vector": 1, "repcap": 2, "helper": 3, "lambda": 4, "host": 5}


def load_annotation_txt(annot_filename):
    """Parse the annotation.txt file into a dictionary.

    Example::

        NAME=chr1;TYPE=host;
        NAME=chr2;TYPE=host;
        NAME=myVector;TYPE=vector;REGION=1795-6553;
        NAME=myCapRep;TYPE=repcap;REGION=1895-5987;
        NAME=myHelper;TYPE=helper;

    :param annot_filename:
        Annotation file following the format indicated above. Only "vector" is required.
        Others optional.
    :return:
        Nested dict like::

        { "myVector": { "label": "vector", "region": (1795, 6553) },
          "myCapRep": { "label": "repcap", "region": (1895, 5987) },
          "myHelper": { "label": "helper", "region": None },
          "chr1":     { "label": "host",   "region": None },
          "chr2":     { "label": "host",   "region": None },
        }

    """
    result = {}
    with open(annot_filename) as inf:
        for line in inf:
            stuff = line.strip()
            m = RX_ANNOT.match(stuff)
            if m is None:
                raise RuntimeError(
                    f"{stuff} is not a valid annotation line! Should follow format "
                    "`NAME=xxxx;TYPE=xxxx;REGION=xxxx;` or `NAME=xxxx;TYPE=xxxx;`. Abort!"
                )

            seq_name = m.group(1)
            ref_label = m.group(2)
            coord_region = (
                None
                if m.group(3) is None
                else tuple(map(int, m.group(3).split("=")[1].split("-")))
            )
            if ref_label in result:
                raise RuntimeError(
                    f"Annotation file has multiple {ref_label} types. Abort!"
                )
            if ref_label not in ANNOT_TYPE_PRIORITIES:
                logging.info(
                    "Nonstandard reference label %s; the known labels are: %s",
                    ref_label,
                    ", ".join(ANNOT_TYPE_PRIORITIES.keys()),
                )
            result[seq_name] = {"label": ref_label, "region": coord_region}
    return result


def anno_dict_to_records(anno_dc: dict):
    """Restructure the dict as BED rows."""
    for seq_name, info in anno_dc.items():
        region = info.get("region")
        if region is not None:
            start1, end = region
            yield (seq_name, start1 - 1, end, info["label"])


if __name__ == "__main__":
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("annotation_txt")
    AP.add_argument("output_bed")
    args = AP.parse_args()

    anno = load_annotation_txt(args.annotation_txt)
    out_rows = anno_dict_to_records(anno)
    with open(args.output_bed, "w") as outf:
        cw = csv.writer(outf, dialect="excel-tab")
        cw.writerows(out_rows)
