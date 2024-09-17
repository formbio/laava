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
    },
    "ss": {
        "alignments": 5489,
        "per_read": 3517,
        "flipflop": 2596,
        "metadata": 1,
        "reference_names": 8,
        "nonmatch": 78330,
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

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_alignments(self, build_dir, name_prefix):
        path = build_dir / f"{name_prefix}.alignments.tsv.gz"
        expected_count = EXPECTED_ROW_COUNTS[name_prefix]["alignments"]
        self.check_row_count(path, expected_count)

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_per_read(self, build_dir, name_prefix):
        tsv_path = build_dir / f"{name_prefix}.per_read.tsv.gz"
        bam_path = Path(f"samples/{name_prefix}.subsample005.bam")
        expected_count = EXPECTED_ROW_COUNTS[name_prefix]["per_read"]

        df = self.check_row_count(tsv_path, expected_count)

        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            primary_count = sum(
                1
                for read in bam_file
                if not read.is_secondary and not read.is_supplementary
            )

        assert (
            len(df) == expected_count == primary_count
        ), f"Expected {expected_count} primary reads, got {len(df)} in TSV and {primary_count} in SAM for {tsv_path}"

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_metadata(self, build_dir, name_prefix):
        path = build_dir / f"{name_prefix}.metadata.tsv"
        expected_count = EXPECTED_ROW_COUNTS[name_prefix]["metadata"]
        self.check_row_count(path, expected_count)

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_reference_names(self, build_dir, name_prefix):
        tsv_path = build_dir / f"{name_prefix}.reference_names.tsv"
        bam_path = build_dir / f"{name_prefix}.tagged.bam"
        expected_count = EXPECTED_ROW_COUNTS[name_prefix]["reference_names"]

        df = self.check_row_count(tsv_path, expected_count)

        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            sn_count = len(bam_file.references)

        assert (
            len(df) == expected_count == sn_count
        ), f"Expected {expected_count} references, got {len(df)} in TSV and {sn_count} in BAM for {tsv_path}"

    @pytest.mark.parametrize("name_prefix", EXPECTED_ROW_COUNTS.keys())
    def test_nonmatch(self, build_dir, name_prefix):
        path = build_dir / f"{name_prefix}.nonmatch.tsv.gz"
        expected_count = EXPECTED_ROW_COUNTS[name_prefix]["nonmatch"]
        self.check_row_count(path, expected_count)

    def test_ss_flipflop(self, build_dir):
        path = build_dir / "ss.flipflop.tsv.gz"
        expected_count = EXPECTED_ROW_COUNTS["ss"]["flipflop"]
        self.check_row_count(path, expected_count)
