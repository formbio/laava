#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { map_reads; make_report; hostgenect; bamqc } from './modules/local/laava'

NO_FILE = file("$projectDir/bin/NO_FILE")

workflow laava {
    take:
    reads
    vector_fa
    packaging_fa
    host_fa
    host_exons
    repcap_name
    vector_bed
    vector_type
    target_gap_threshold
    max_allowed_outside_vector
    max_allowed_missing_flanking
    flipflop_name
    flipflop_fa

    main:
    map_reads(
        reads
        .combine(vector_fa)
        .combine(packaging_fa)
        .combine(host_fa)
        .combine(repcap_name))
    make_report(
        map_reads.out.mapped_reads
        .combine(vector_bed)
        .combine(vector_type)
        .combine(target_gap_threshold)
        .combine(max_allowed_outside_vector)
        .combine(max_allowed_missing_flanking)
        .combine(flipflop_name)
        .combine(flipflop_fa))
    if (host_exons) {
      hostgenect(map_reads.out.bam.combine(host_exons))
      hostct=hostgenect.out
    } else {
      hostct=Channel.of(NO_FILE)
    }
    bamqc(map_reads.out.bam)
    emit:
    sam = map_reads.out.mapped_reads
    per_read_tsv = make_report.out.per_read_tsv
    summary_tsv = make_report.out.summary_tsv
    nonmatch_stat_tsvgz = make_report.out.nonmatch_stat_tsvgz
    tagged_bam = make_report.out.tagged_bam
    subtype_bams = make_report.out.subtype_bams
    subtype_bais = make_report.out.subtype_bais
    flipflop_assignments_tsv = make_report.out.flipflop_assignments_tsv
    flipflop_bams = make_report.out.flipflop_bams
    alignments_tsv = make_report.out.alignments_tsv
    readsummary_tsv = make_report.out.readsummary_tsv
    sequence_error_tsv = make_report.out.sequence_error_tsv
    flipflop_tsv = make_report.out.flipflop_tsv
    rdata = make_report.out.rdata
    qc=bamqc.out.qc
    hostct=hostct
}

workflow {
    if (params.inputdir) {
      seqfiles = Channel.fromPath(["${params.inputdir}/*.fastq.gz","${params.inputdir}/*.fq.gz","${params.inputdir}/*.fastq","${params.inputdir}/*.fq"])
        .map { file -> tuple( file.getName().split(/\.fq|\.fastq/)[0],file ) }
    } else {
      seqfiles = Channel.fromPath(params.seq_reads)
        .map { file -> tuple(file.getName().split(/\.fq|\.fastq/)[0], file) }
    }
    laava(
        seqfiles,
        Channel.fromPath(params.vector_fa),
        params.packaging_fa ? Channel.fromPath(params.packaging_fa) : Channel.of(NO_FILE),
        params.host_fa ? Channel.fromPath(params.host_fa) : Channel.of(NO_FILE),
        params.host_exons ? Channel.fromPath(params.host_exons) : false,
        Channel.of(params.repcap_name),
        Channel.fromPath(params.vector_bed),
        Channel.of(params.vector_type),
        Channel.of(params.target_gap_threshold),
        Channel.of(params.max_allowed_outside_vector),
        Channel.of(params.max_allowed_missing_flanking),
        Channel.of(params.flipflop_name),
        params.flipflop_fa ? Channel.fromPath(params.flipflop_fa) : Channel.of(NO_FILE)
    )
}
