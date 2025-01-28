## Summary

This workflow can be used with PacBio long-read sequencing to characterize
adeno-associated virus (AAV) products by examining alignment orientation and
coverage over annotated sequence regions. It produces a report summarizing the
analysis results and supporting tabular data.

## Workflow Diagram

``` mermaid
graph TD

seqreads(Sequencing reads) --> map_reads(Read alignment)
refseqs(Reference sequences) --> map_reads
map_reads --> classify(Read classification)
anno(Annotation) --> classify
meta(Sample metadata) --> classify
classify --> aggtables(Summary tables)
classify --> report(Report)
aggtables --> report(Report)
```

## Methods

This workflow performs an automated analysis of long-read sequencing data from
adeno-associated virus (AAV) products. The sequencing data should be from a
PacBio sequencer (Sequel II, Sequel IIe or Revio) run with on-instrument AAV
mode enabled, or equivalent circular consensus sequencing (CCS) reads (Travers
et al., 2010)

In this analysis, reads are aligned to the given AAV, packaging, and host reference
sequences using Minimap2 (Li, 2018).
The reference sequences for each primary alignment and its orientation are counted and
summarized to assign read type classifications, including vector, non-vector, and
chimeric reads.
For reads assigned to the AAV vector, the primary alignment coordinates are compared to
the annotated vector region in the reference sequence, which comprises the left and
right ITRs and the genomic region between them, to assign each read to a subtype
classification.
Sequence variants relative to the vector reference sequence are determined directly from
each read's alignment, specifically the CIGAR string indicating insertions, deletions,
mismatches, and gaps.

Finally, a report is generated with relevant quality metrics and analysis
results in both HTML and PDF formats.
