## Inputs

### Files

- Sequencing reads (unmapped BAM)
  - PacBio AAV sequencing read set, as HiFi/CCS reads in unaligned BAM format. The
    PacBio instrument should be run in AAV mode.
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


### Options

- AAV Serotype
  - Serotype name for ITR flip/flop analysis, e.g. "AAV2", or leave blank to skip.
- Thresholds for "full" read classification
    - Target gap threshold
    - Max allowed outside vector
    - Max allowed missing flanking
