#!/usr/bin/env python3
"""LAAVA Runner Script - Unified command line interface for LAAVA analysis pipeline."""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

# Add src directory to Python path
src_dir = str(Path(__file__).parent.parent)
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

# Import the modules
import importlib.util  # noqa: E402

get_reference_names = importlib.import_module("src.get_reference_names")
guess_vector_type_length = importlib.import_module("src.guess_vector_type_length")
summarize_alignment = importlib.import_module("src.summarize_alignment")

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="LAAVA Analysis Pipeline Runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        "--input-bam", required=True, help="Input BAM file containing sequencing reads"
    )
    parser.add_argument(
        "--vector-fasta",
        required=True,
        help="Vector reference sequence in FASTA format",
    )
    parser.add_argument(
        "--annotation-bed", required=True, help="Vector annotation in BED format"
    )
    parser.add_argument(
        "--itr-label", required=True, help="Label used for ITR in annotation"
    )

    # Optional arguments
    parser.add_argument(
        "--packaging-fasta", help="Packaging sequence reference in FASTA format"
    )
    parser.add_argument("--host-fasta", help="Host genome reference in FASTA format")
    parser.add_argument(
        "--output-dir",
        default="./laava_results",
        help="Output directory for analysis results",
    )
    parser.add_argument("--sample-id", default="sample", help="Sample identifier")
    parser.add_argument("--sample-name", help="Human-readable sample name")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    return parser.parse_args()


def setup_output_directory(output_dir: str) -> Path:
    """Create and validate output directory."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    return out_path


def run_pipeline(args: argparse.Namespace) -> int:
    """Execute the LAAVA analysis pipeline."""
    try:
        # Setup output directory
        output_dir = setup_output_directory(args.output_dir)

        # Step 1: Analyze vector type
        logger.info("Analyzing vector type...")
        cassette_len = guess_vector_type_length.length_from_annotation(
            args.annotation_bed, [args.itr_label]
        )
        logger.info(f"Expression cassette length is {cassette_len}")
        vector_type = (
            "sc" if cassette_len <= guess_vector_type_length.SC_MAX_THRESHOLD else "ss"
        )
        logger.info(f"Determined vector type: {vector_type}")

        # Step 2: Generate reference names
        logger.info("Generating reference names...")
        ref_names_file = output_dir / f"{args.sample_id}.reference_names.tsv"
        ref_names_args = argparse.Namespace(
            vector=args.vector_fasta,
            packaging=args.packaging_fasta,
            host=args.host_fasta,
            output=str(ref_names_file),
            repcap_name=None,
            helper_name=None,
            lambda_name=None,
        )
        get_reference_names._main(ref_names_args)

        # Step 3: Map reads
        logger.info("Mapping reads...")

        # Set up environment for the shell script
        env = os.environ.copy()
        env["PATH"] = (
            f"{Path(__file__).parent}:{env['PATH']}"  # Add src dir to PATH for Python scripts
        )

        # Get number of CPUs for parallel processing
        try:
            import multiprocessing

            threads = str(multiprocessing.cpu_count())
        except (NotImplementedError, ValueError):
            threads = "1"
        env["THREADS"] = threads

        # Prepare command
        map_script = Path(__file__).parent / "map_reads.sh"
        if not map_script.exists():
            raise FileNotFoundError(f"Map reads script not found: {map_script}")

        # Build command with proper working directory
        cmd = [
            "bash",
            str(map_script),
            args.sample_id,
            str(Path(args.input_bam).absolute()),
            str(Path(args.vector_fasta).absolute()),
        ]
        if args.packaging_fasta:
            cmd.append(str(Path(args.packaging_fasta).absolute()))
        if args.host_fasta:
            cmd.append(str(Path(args.host_fasta).absolute()))

        # Run mapping command in output directory
        logger.debug(f"Running command: {' '.join(cmd)}")
        try:
            subprocess.run(
                cmd,
                check=True,
                cwd=output_dir,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logger.error("Map reads command failed:")
            logger.error(f"stdout:\n{e.stdout}")
            logger.error(f"stderr:\n{e.stderr}")
            raise

        # Verify output file exists
        mapped_reads = output_dir / f"{args.sample_id}.sort_by_name.sam"
        if not mapped_reads.exists():
            raise FileNotFoundError(
                f"Expected mapped reads file not found: {mapped_reads}"
            )

        # Step 4: Analyze alignments and generate report
        logger.info("Analyzing alignments and generating report...")
        analysis_args = argparse.Namespace(
            sam_filename=str(mapped_reads),
            annotation_bed=args.annotation_bed,
            reference_names=str(ref_names_file),
            itr_labels=[args.itr_label],
            output_prefix=str(output_dir / args.sample_id),
            sample_id=args.sample_id,
            vector_type=vector_type,
            cpus=1,
            debug=args.verbose,
        )
        summarize_alignment.main(analysis_args)

        logger.info(f"Analysis complete. Results available in: {output_dir}")
        return 0

    except Exception as e:
        logger.error(f"Pipeline failed: {e!s}")
        if args.verbose:
            logger.exception("Detailed error traceback:")
        return 1


def main():
    """Main entry point for the LAAVA runner script."""
    args = parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    sys.exit(run_pipeline(args))


if __name__ == "__main__":
    main()
