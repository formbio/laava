"""Test the outputs of the local LAAVA commands.

Usage: pytest test_outputs.py

Before running this test suite, first run the LAAVA commands using the accompanying
Makefile. That will generate the outputs this test suite looks for.
"""

from pathlib import Path

import pandas as pd
import pysam
import pytest

EXPECTED_ROW_COUNTS = {
    "sc": {
        "alignments": 3349,
        "per_read": 1687,
        "metadata": 1,
        "reference_names": 8,
        "nonmatch": 52620,
        "agg_ref_type": 4,
        "agg_subtype": 9,
    },
    "ss": {
        "alignments": 5489,
        "per_read": 3517,
        "flipflop": 2576,
        "metadata": 1,
        "reference_names": 8,
        "nonmatch": 78330,
        "agg_ref_type": 9,
        "agg_subtype": 13,
        "agg_flipflop": 15,
    },
}

BUILD_DIR = "build"


class TestCompareTSVs:
    @pytest.fixture(scope="class")
    def build_dir(self):
        return Path(BUILD_DIR)

    def check_row_count(self, path, expected_count):
        df = pd.read_csv(path, sep="\t")
        assert (
            len(df) == expected_count
        ), f"Expected {expected_count} rows, got {len(df)} for {path}"
        return df

    def check_tsv_row_count(self, build_dir, name_prefix, key, suffix):
        tsv_path = build_dir / f"{name_prefix}.{key}.{suffix}"
        expected_count = EXPECTED_ROW_COUNTS[name_prefix][key]
        df = self.check_row_count(tsv_path, expected_count)
        return tsv_path, expected_count, len(df)

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_metadata(self, build_dir, name_prefix):
        self.check_tsv_row_count(build_dir, name_prefix, "metadata", "tsv")

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_alignments(self, build_dir, name_prefix):
        self.check_tsv_row_count(build_dir, name_prefix, "alignments", "tsv.gz")

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_nonmatch(self, build_dir, name_prefix):
        self.check_tsv_row_count(build_dir, name_prefix, "nonmatch", "tsv.gz")

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_agg_ref_type(self, build_dir, name_prefix):
        self.check_tsv_row_count(build_dir, name_prefix, "agg_ref_type", "tsv")

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_agg_subtype(self, build_dir, name_prefix):
        self.check_tsv_row_count(build_dir, name_prefix, "agg_subtype", "tsv")

    # Only ssAAV
    def test_ss_flipflop(self, build_dir):
        self.check_tsv_row_count(build_dir, "ss", "flipflop", "tsv.gz")

    def test_ss_agg_flipflop(self, build_dir):
        self.check_tsv_row_count(build_dir, "ss", "agg_flipflop", "tsv")

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_per_read(self, build_dir, name_prefix):
        tsv_path, expected_count, row_count = self.check_tsv_row_count(
            build_dir, name_prefix, "per_read", "tsv.gz"
        )
        bam_path = Path(f"samples/{name_prefix}.subsample005.bam")
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            primary_count = sum(
                1
                for read in bam_file
                if not read.is_secondary and not read.is_supplementary
            )
        assert (
            row_count == expected_count == primary_count
        ), f"Expected {expected_count} primary reads, got {row_count} in TSV and {primary_count} in SAM for {tsv_path}"

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_reference_names(self, build_dir, name_prefix):
        tsv_path, expected_count, row_count = self.check_tsv_row_count(
            build_dir, name_prefix, "reference_names", "tsv"
        )
        bam_path = build_dir / f"{name_prefix}.tagged.bam"
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            sn_count = len(bam_file.references)
        assert (
            row_count == expected_count == sn_count
        ), f"Expected {expected_count} references, got {row_count} in TSV and {sn_count} in BAM for {tsv_path}"
