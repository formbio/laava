"""Golden master verification tests for LAAVA pipeline outputs.

This is a "naive catch-all" test suite that detects ANY content changes in pipeline
outputs by comparing file content hashes against stored golden masters.

PURPOSE:
========
- Catch subtle algorithm changes that don't affect row counts
- Detect rounding/precision changes in statistics  
- Flag any modification to data values or column ordering
- Force developers to make conscious decisions about output changes

WORKFLOW:
=========
1. Pipeline generates outputs via Makefile
2. Tests hash the content of critical output files
3. Compare against stored golden hashes
4. FAIL if ANY content has changed
5. Developer must either:
   - Fix the bug that caused the change, OR
   - Update the golden hashes if change was intentional

This complements other test types:
- test_outputs.py: Structural validation (row counts, formats)
- test_integration_duplicate_removal.py: Business logic validation
- test_unit_laava_modules.py: Direct code testing
- test_output_gold_verification.py: Comprehensive content verification
"""

import hashlib
import gzip
import pytest
from pathlib import Path
import pandas as pd


class TestGoldVerification:
    """Golden master tests that detect ANY content changes in pipeline outputs.
    
    These tests are intentionally strict - they will fail if ANYTHING changes
    in the output files, forcing developers to consciously decide whether
    changes are bugs (fix the code) or features (update the hashes).
    """

    # GOLDEN HASHES: Update these when changes are intentional
    # To update: run the pipeline, then copy the new hashes from test failures
    GOLDEN_HASHES = {
        # SC (Self-Complementary) pipeline outputs
        "sc": {
            "nonmatch": "cf9cf24d7efdc43086cf47e6aa0977b5bfbd22b2e7915b59b79d291791a5dd4d",
            "alignments": "28e41d7f23699fe4e0a7dfbde2039c9f7de8bf249d2c771f98c3cdac42c8c78f", 
            "per_read": "fd694fb03859deb296f95fe4886ba04642d5a60e689fbb87cbddd390544d076a",
            "agg_ref_type": "7a7635f8896cecd142438f6ceb868939ca642593f2df2b6dfb19f58c38b28d37",
            "agg_subtype": "734b36a42ad47c00f4e61c0194ba129b74a74b9a10d507427f3a6f600e0bd87e",
            "metadata": "0f9d78923ac48bec7b204900d45d57be62f2d0c241fb5a769f3691f8d61a3dcc",
        },
        # SS (Single-Stranded) pipeline outputs  
        "ss": {
            "nonmatch": "7547a6278e470637ddd14c40b6b4c72c27983966af8c821333237d3868a644e5",
            "alignments": "9d57fa00c1f60b1bf84281ce21bf8e8deed9b4a0e50e0bfdb37984c2a9558b0e",
            "per_read": "26aefedb8a2d0b747015d999bff9e35e5acd7f040f5f89efeb40289598ff666a", 
            "agg_ref_type": "1fe190033bed259a61f79e6bad8b01151132e6fa58ddc99b50fd5db2e1bb48f1",
            "agg_subtype": "c58e2d0f04e5e3ac0a8a5eb94e5e885719bca6d4a019f288b8b0f0089b715552",
            "agg_flipflop": "0fbf0b00fc757989160de05e238032a07fc5108de3de59fd62e6d9267fbb8521",
            "flipflop": "6b80dbe7a6b26c0b66944ff6dd08947cde8c1d142553424e51f0e85d7fb96704",
            "metadata": "374c65f7a71d1607263d54fdfcee879de390cf1e718a5a6c3c1a71d3e32e3729",
        },
        # TC-GIA-012 pipeline outputs
        "tc-gia-012": {
            "nonmatch": "d6473fd0260da7858a155a8d00756a87553353a4c47f836f4e57f2c4f4afb673",
            "alignments": "682fe02480cc81a7725b26936bc1b38b10370506501e97f302e2d73a0defcef6",
            "per_read": "828211a905e30653312d907118258ac2ced59b986931eb036f04fa1c35ec3890",
            "agg_ref_type": "2149128cd3454adc1ad2a00633246984df142127e0dc58a739bdf27c5e6b934c",
            "agg_subtype": "eb80b89e4d7357fda99cd802544c39ee9f8b310c1026cd30983faa6bfa0903fb",
            "agg_flipflop": "ae361932caddc73c33972f80341d0906130e7386b915da4435fb292e22f4a4a3",
            "flipflop": "de19ecfca2f2167dc672f5e9a608469a6515f205394d99c53ae257465f6b7aae",
            "metadata": "acbd7bef18c335d4c5ca0ffa3aae89c6cb3b04aac07ab5214af8898c70572304",
        }
    }

    @pytest.fixture(scope="class")
    def build_dir(self):
        """Get the build directory path."""
        return Path("build")

    def hash_file_content(self, file_path):
        """Calculate SHA256 hash of file content.
        
        For gzipped files, decompresses and hashes the actual content to avoid
        timestamp differences in gzip headers that cause false positives.
        For TSV files, normalizes the content to handle minor formatting differences.
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Output file not found: {file_path}")
        
        # For gzipped files, decompress and hash the actual content
        if file_path.suffix == '.gz':
            with gzip.open(file_path, 'rt') as f:
                content = f.read().encode('utf-8')
        else:
            # For TSV files, normalize content to handle formatting differences
            try:
                # Read as DataFrame and convert back to standardized TSV
                df = pd.read_csv(file_path, sep='\t')
                content = df.to_csv(sep='\t', index=False).encode('utf-8')
            except Exception:
                # Fallback to raw content if pandas can't read it
                with open(file_path, 'rb') as f:
                    content = f.read()
        
        return hashlib.sha256(content).hexdigest()

    def verify_file_hash(self, build_dir, sample_type, file_key, file_suffix):
        """Verify a single file's content hash against golden master.
        
        Args:
            build_dir: Path to build directory
            sample_type: 'sc' or 'ss'
            file_key: Key in GOLDEN_HASHES dict
            file_suffix: File extension (e.g., 'tsv', 'tsv.gz')
        """
        file_path = build_dir / f"{sample_type}.{file_key}.{file_suffix}"
        
        # Calculate current hash
        current_hash = self.hash_file_content(file_path)
        
        # Get expected hash
        expected_hash = self.GOLDEN_HASHES[sample_type][file_key]
        
        # Check if this is a placeholder (first run)
        if expected_hash.startswith("PLACEHOLDER_"):
            pytest.skip(f"Golden hash not set for {sample_type}.{file_key} - "
                       f"current hash: {current_hash}")
        
        # Verify hash matches
        assert current_hash == expected_hash, (
            f"Content changed in {file_path}!\n"
            f"Expected hash: {expected_hash}\n"
            f"Current hash:  {current_hash}\n\n"
            f"If this change was intentional, update GOLDEN_HASHES['{sample_type}']['{file_key}'] "
            f"to '{current_hash}'"
        )

    # SC (Self-Complementary) Pipeline Tests
    def test_sc_nonmatch_content_unchanged(self, build_dir):
        """Verify SC nonmatch.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "sc", "nonmatch", "tsv.gz")

    def test_sc_alignments_content_unchanged(self, build_dir):
        """Verify SC alignments.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "sc", "alignments", "tsv.gz")

    def test_sc_per_read_content_unchanged(self, build_dir):
        """Verify SC per_read.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "sc", "per_read", "tsv.gz")

    def test_sc_agg_ref_type_content_unchanged(self, build_dir):
        """Verify SC agg_ref_type.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "sc", "agg_ref_type", "tsv")

    def test_sc_agg_subtype_content_unchanged(self, build_dir):
        """Verify SC agg_subtype.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "sc", "agg_subtype", "tsv")

    def test_sc_metadata_content_unchanged(self, build_dir):
        """Verify SC metadata.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "sc", "metadata", "tsv")

    # SS (Single-Stranded) Pipeline Tests
    def test_ss_nonmatch_content_unchanged(self, build_dir):
        """Verify SS nonmatch.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "nonmatch", "tsv.gz")

    def test_ss_alignments_content_unchanged(self, build_dir):
        """Verify SS alignments.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "alignments", "tsv.gz")

    def test_ss_per_read_content_unchanged(self, build_dir):
        """Verify SS per_read.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "per_read", "tsv.gz")

    def test_ss_agg_ref_type_content_unchanged(self, build_dir):
        """Verify SS agg_ref_type.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "agg_ref_type", "tsv")

    def test_ss_agg_subtype_content_unchanged(self, build_dir):
        """Verify SS agg_subtype.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "agg_subtype", "tsv")

    def test_ss_agg_flipflop_content_unchanged(self, build_dir):
        """Verify SS agg_flipflop.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "agg_flipflop", "tsv")

    def test_ss_flipflop_content_unchanged(self, build_dir):
        """Verify SS flipflop.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "flipflop", "tsv.gz")

    def test_ss_metadata_content_unchanged(self, build_dir):
        """Verify SS metadata.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "ss", "metadata", "tsv")

    # TC-GIA-012 Pipeline Tests
    def test_tc_gia_012_nonmatch_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 nonmatch.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "nonmatch", "tsv.gz")

    def test_tc_gia_012_alignments_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 alignments.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "alignments", "tsv.gz")

    def test_tc_gia_012_per_read_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 per_read.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "per_read", "tsv.gz")

    def test_tc_gia_012_agg_ref_type_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 agg_ref_type.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "agg_ref_type", "tsv")

    def test_tc_gia_012_agg_subtype_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 agg_subtype.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "agg_subtype", "tsv")

    def test_tc_gia_012_agg_flipflop_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 agg_flipflop.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "agg_flipflop", "tsv")

    def test_tc_gia_012_flipflop_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 flipflop.tsv.gz content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "flipflop", "tsv.gz")

    def test_tc_gia_012_metadata_content_unchanged(self, build_dir):
        """Verify TC-GIA-012 metadata.tsv content is unchanged."""
        self.verify_file_hash(build_dir, "tc-gia-012", "metadata", "tsv")

    # Utility Methods
    def print_current_hashes(self, build_dir):
        """Helper method to print current hashes for updating GOLDEN_HASHES.
        
        Run this manually when you want to update the golden masters:
        pytest test_output_gold_verification.py::TestGoldVerification::print_current_hashes -s
        """
        print("\n" + "="*60)
        print("CURRENT FILE HASHES (for updating GOLDEN_HASHES)")
        print("="*60)
        
        for sample_type in ["sc", "ss"]:
            print(f"\n'{sample_type}': {{")
            
            file_specs = [
                ("nonmatch", "tsv.gz"),
                ("alignments", "tsv.gz"), 
                ("per_read", "tsv.gz"),
                ("agg_ref_type", "tsv"),
                ("agg_subtype", "tsv"),
                ("metadata", "tsv"),
            ]
            
            # Add SS-specific files
            if sample_type == "ss":
                file_specs.extend([
                    ("agg_flipflop", "tsv"),
                    ("flipflop", "tsv.gz"),
                ])
            
            for file_key, file_suffix in file_specs:
                file_path = build_dir / f"{sample_type}.{file_key}.{file_suffix}"
                if file_path.exists():
                    current_hash = self.hash_file_content(file_path)
                    print(f"    '{file_key}': '{current_hash}',")
                else:
                    print(f"    '{file_key}': 'FILE_NOT_FOUND',")
            
            print("},")
        
        print("\n" + "="*60)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
