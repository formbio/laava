"""Test the additional (additive) outputs: appended columns, the corrected map_iden and all_alignments.tsv.gz.

Usage: pytest test_additive_outputs.py (after running the Makefile, like test_outputs.py)
"""

from pathlib import Path

import pandas as pd
import pysam
import pytest

SAMPLES = ["sc", "ss"]
BUILD = Path("build")


def load(name):
    return pd.read_csv(BUILD / name, sep="\t", dtype=str, keep_default_na=False)


@pytest.mark.parametrize("sample", SAMPLES)
def test_appended_columns(sample):
    aln = load(f"{sample}.alignments.tsv.gz")
    assert list(aln.columns[-4:]) == ["mapq", "s1", "s2", "de"]
    assert aln.columns[10] == "map_iden"  # original columns keep their positions
    per_read = load(f"{sample}.per_read.tsv.gz")
    assert list(per_read.columns[-2:]) == ["n_alignments", "labels_involved"]


@pytest.mark.parametrize("sample", SAMPLES)
def test_map_iden_is_identity_over_alignment_block(sample):
    """map_iden = 1 - (X+I+D+N) / (=+X+I+D+N) per alignment (the minimap2 PAF matches/block definition)."""
    aln = load(f"{sample}.alignments.tsv.gz")
    aln = aln[aln["is_mapped"] == "Y"].copy()
    aln["key"] = aln.read_id + "|" + aln.map_start0 + "|" + aln.map_end1 + "|" + aln.is_supp + "|" + aln.map_strand
    expected = {}
    for r in pysam.AlignmentFile(str(BUILD / f"{sample}.tagged.bam"), "rb", check_sq=False):
        if r.is_unmapped:
            continue
        c = {k: 0 for k in "=XIDN"}
        for op, n in r.cigartuples:
            ch = "MIDNSHP=X"[op]
            if ch in c:
                c[ch] += n
        block = sum(c.values())
        key = f"{r.qname}|{r.reference_start}|{r.reference_end}|{'Y' if r.is_supplementary else 'N'}|{'-' if r.is_reverse else '+'}"
        expected[key] = 1 - (c["X"] + c["I"] + c["D"] + c["N"]) / block
    got = aln.set_index("key")["map_iden"].astype(float)
    assert set(got.index) <= set(expected)
    assert max(abs(got[k] - expected[k]) for k in got.index) < 1e-9


@pytest.mark.parametrize("sample", SAMPLES)
def test_all_alignments_file(sample):
    aln = load(f"{sample}.alignments.tsv.gz")
    allal = load(f"{sample}.all_alignments.tsv.gz")
    assert list(allal.columns) == list(aln.columns) + ["kept"]
    kept = allal[allal["kept"] == "Y"].drop(columns="kept").reset_index(drop=True)
    assert kept.equals(aln.reset_index(drop=True))
    assert set(allal["kept"]) <= {"Y", "N"}


@pytest.mark.parametrize("sample", SAMPLES)
def test_n_alignments_matches_all_alignments(sample):
    per_read = load(f"{sample}.per_read.tsv.gz").set_index("read_id")
    allal = load(f"{sample}.all_alignments.tsv.gz")
    mapped = allal[allal["is_mapped"] == "Y"].groupby("read_id").size()
    n = per_read["n_alignments"].astype(int)
    assert n[n > 0].sort_index().equals(mapped.sort_index())
    assert (n[per_read["has_supp"] == "Y"] >= 2).all()
    unmapped = per_read[per_read["has_primary"] != "Y"]
    assert (unmapped["labels_involved"] == "NA").all()
