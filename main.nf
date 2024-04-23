#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { map_reads; make_report } from './modules/local/laava'

NO_FILE = file("$projectDir/bin/NO_FILE")

workflow {
    result = map_reads(
        file(params.seq_reads).getSimpleName(),
        Channel.fromPath(params.seq_reads),
        Channel.fromPath(params.vector_fa),
        params.packaging_fa ? Channel.fromPath(params.packaging_fa) : NO_FILE,
        params.host_fa ? Channel.fromPath(params.host_fa): NO_FILE,
        params.repcap_name,
    )
    make_report(
        Channel.fromPath(params.vector_bed),
        result.reference_names,
        result.sample_name,
        params.target_gap_threshold,
        params.max_allowed_outside_vector,
        params.max_allowed_missing_flanking,
        result.mapped_reads,
        params.flipflop_name,
        params.flipflop_fa ? Channel.fromPath(params.flipflop_fa): NO_FILE,
    )
}
