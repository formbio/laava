"""`map_reads.sh` keeps the input BAM's header as `<sample>.source_header.txt` (BAM input only), before `samtools fastq` drops it.

Needs `minimap2`, `samtools` and this repository's `src/` scripts on the PATH (like the rest of the suite inside the LAAVA image). Runs on a tiny synthetic reference and reads.
"""

import os
import random
import shutil
import subprocess
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None or shutil.which("samtools") is None, reason="minimap2 and samtools are required")

HEADER = (
    "@HD\tVN:1.6\tSO:unknown\n"
    "@RG\tID:abc123/0--0\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;BASECALLERVERSION=5.0\tPU:m64248e_241116_225550\tSM:sample\n"
    "@PG\tID:ccs\tPN:ccs\tVN:8.0.0\n"
)


def _make_inputs(tmp: Path):
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(3000))
    (tmp / "vector.fa").write_text(">vec\n" + ref + "\n")
    reads = [ref[i * 100 : i * 100 + 800] for i in range(8)]
    sam = HEADER + "".join(f"m64248e_241116_225550/{i}/ccs\t4\t*\t0\t0\t*\t*\t0\t0\t{r}\t{'I' * len(r)}\tRG:Z:abc123/0--0\n" for i, r in enumerate(reads))
    (tmp / "in.sam").write_text(sam)
    subprocess.run(["samtools", "view", "-b", "-o", str(tmp / "reads.bam"), str(tmp / "in.sam")], check=True)
    (tmp / "reads.fq").write_text("".join(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n" for i, r in enumerate(reads)))


def _run(tmp: Path, reads: str):
    env = dict(os.environ, PATH=f"{SRC}:{os.environ['PATH']}")
    subprocess.run(["bash", str(SRC / "map_reads.sh"), "s1", reads, "vector.fa", "", "", "", "", ""], cwd=tmp, env=env, check=True, capture_output=True)


def test_a_bam_input_keeps_its_header_verbatim(tmp_path):
    _make_inputs(tmp_path)
    _run(tmp_path, "reads.bam")
    kept = (tmp_path / "s1.source_header.txt").read_text()
    assert kept == subprocess.run(["samtools", "view", "-H", "--no-PG", str(tmp_path / "reads.bam")], capture_output=True, text=True, check=True).stdout
    assert "@RG\tID:abc123/0--0" in kept and "BINDINGKIT=102-739-100" in kept
    # the aligned BAM does not have it (this is the loss the file works around)
    assert "@RG" not in subprocess.run(["samtools", "view", "-H", str(tmp_path / "s1.sort_by_pos.bam")], capture_output=True, text=True, check=True).stdout


def test_a_fastq_input_writes_no_header_file_and_still_runs(tmp_path):
    _make_inputs(tmp_path)
    _run(tmp_path, "reads.fq")
    assert not (tmp_path / "s1.source_header.txt").exists()
    assert (tmp_path / "s1.sort_by_pos.bam").is_file()
