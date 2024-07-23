## Inputs

### Files

For a single sample:

- Sequencing reads file (unmapped BAM or FASTQ)
  - PacBio AAV sequencing read set, as HiFi/CCS reads in unaligned BAM format. The
    PacBio instrument should be run in AAV mode.
- Sample unique ID
  - A unique identifier to distinguish the sample from others in the same sequencing run. This may be the barcode number used in the sample sheet, or a generated UUID.",
- Sample display name
  - Sample name to display in the output report. Not required to be unique; e.g.
    replicates may share a display name but have different unique IDs.

For a folder of multiple samples (post-demux):

- Sequencing reads folder (folder of unmapped BAM or FASTQ)
    - As above.
- Sample metadata file (TSV)
    - Unique ID and display name for each sample in a two-column, tab-separated file.
      The column headers must be 'sample_unique_id' (e.g. barcode number, unique within
      a sequencing run) and 'sample_display' name (to display in the report, not
      required to be unique).

The remaining parameters are shared by all samples:

- Vector annotation (BED)
  - Annotated vector construct region coordinates in 4-column BED format. Must include 1
    region with the label 'vector', indicating the payload region that spans ITR to ITR.
    May also include additional labeled regions, e.g.  for ITR, promoter, and CDS
    regions, but these will not affect the output.
- Reference sequences (FASTA):
  - Vector plasmid, as a single-record FASTA.
  - Packaging sequences (optional) -- helper and rep/cap plasmids, and other sequences
    of interest (e.g. Lambda), as a multi-record FASTA.
  - Host genome (optional but recommended), as a multi-record FASTA. Best to include
    only the canonical chromosomes and not the alternative contigs.
- Flip/flop ITR sequences (FASTA; optional)
  - AAV2 sequences are built in and available by default; provide custom sequences here
    to use another serotype (and use a different value for the "AAV Serotype" option).


### Labels

These string inputs are used to guide the handling of the file inputs above.

- ITR labels used in the vector annotation BED file (itr_label_1, itr_label_2)
- Sequence IDs used in the packaging FASTA file (repcap_name, helper_name, lambda_name)
- ITR AAV Serotype name for flip/flop analysis (flipflop_name) -- use to select the
  built-in AAV2 flip/flop sequences.

### Options

- Vector type
  - Whether the vector construct is designed to be single-stranded AAV (ssAAV), self-complementary AAV (scAAV), or unknown/unspecified.
- ITR AAV Serotype
  - Serotype name for ITR flip/flop analysis, e.g. "AAV2", or leave blank to skip.
- Thresholds for "full" read classification
    - Target gap threshold
    - Max allowed outside vector
    - Max allowed missing flanking
