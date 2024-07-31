## Inputs

### Files

For a single sample:

- Sequencing reads file (unmapped BAM or FASTQ)
  - PacBio AAV sequencing read set, as HiFi/CCS reads in unaligned BAM format. The
    PacBio instrument should be run in AAV mode.
- Sample unique ID
  - A unique identifier to distinguish the sample from others in the same sequencing
    run. This may be the barcode number used in the sample sheet, or a generated UUID.",
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
  - Annotated vector construct region coordinates in 4-column UCSC BED format.
  - This file must indicate the transgene/payload region via either two labeled Inverted
    Terminal Repeat (ITR) regions (see itr_label_1 and itr_label_2 below) or, as a
    legacy mode, one region with the label 'vector', spanning both ITRs (inclusive).
  - May also include additional labeled regions, e.g. for promoter and CDS regions;
    these will be ignored and will not affect the output.
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

- ITR labels used in the vector annotation BED file: itr_label_1, itr_label_2
    - These are case-sensitive and must match exactly.
    - The order does not matter: LAAVA will check for the presence of both labels in the
      annotation BED file, and treat the first match as the left (5') ITR and the second
      as the right (3') ITR. If the annotation uses the same label for both, e.g. 'ITR',
      you only need to specify itr_label_1 as 'ITR', and that label will correctly match
      both ITR regions in the annotation BED.
    - If both fields are left blank, LAAVA will fall back to the legacy mode of looking
      for a region labeled "vector" instead.
- Sequence IDs used in the packaging FASTA file: repcap_name, helper_name, lambda_name
    - These are case-sensitive and must match exactly.
    - If not specified, the reference sequence IDs in the packaging FASTA file will be
      counted and reported with their original names. If specified but not found in the
      packaging FASTA file, LAAVA will raise an error.
    - If repcap_name is specified and found in the packaging FASTA, certain plots
      related to RepCap read alignments will be included in the output report.
- ITR AAV Serotype name (flipflop_name)
    - Selects a set of built-in ITR sequences for flip/flop analysis. Currently, the
      only built-in set of ITR sequences is for the AAV2 serotype.
    - Specifying 'AAV2' here is equivalent to providing the same sequences as a FASTA
      file via the "flip/flop ITR sequences" input field above.


### Options

- Vector type
  - Whether the vector construct is designed to be single-stranded AAV (ssAAV),
    self-complementary AAV (scAAV), or unknown/unspecified.
- ITR AAV Serotype
  - Serotype name for ITR flip/flop analysis, e.g. "AAV2", or leave blank to skip.
- Thresholds for "full" read classification
    - Target gap threshold
    - Max allowed outside vector
    - Max allowed missing flanking
