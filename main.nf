#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { map_reads; make_report } from './modules/local/laava'

NO_FILE = file("$projectDir/bin/NO_FILE")
NO_FILE2 = file("$projectDir/bin/NO_FILE2")

workflow laava {
    take:
    reads
    vector_fa
    packaging_fa
    host_fa
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
        map_reads.out.mapped_sam
        .combine(vector_bed)
        .combine(vector_type)
        .combine(target_gap_threshold)
        .combine(max_allowed_outside_vector)
        .combine(max_allowed_missing_flanking)
        .combine(flipflop_name)
        .combine(flipflop_fa))

    emit:
    mapped_sam = map_reads.out.mapped_sam
    mapped_bam = map_reads.out.mapped_bam
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
}

workflow {
    if (params.seq_reads_folder) {
        seq_source = Channel.fromPath([
                "${params.seq_reads_folder}/*.fastq.gz",
                "${params.seq_reads_folder}/*.fq.gz",
                "${params.seq_reads_folder}/*.fastq",
                "${params.seq_reads_folder}/*.fq"])
    } else {
        seq_source = Channel.fromPath(params.seq_reads_file)
    }
    seq_files = seq_source.map {
        file -> tuple(file.getName().split(/\.bam|\.fq|\.fastq/)[0], file)
    }

    laava(
        seq_files,
        Channel.fromPath(params.vector_fa),
        params.packaging_fa ? Channel.fromPath(params.packaging_fa) : Channel.of(NO_FILE),
        params.host_fa ? Channel.fromPath(params.host_fa) : Channel.of(NO_FILE2),
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
