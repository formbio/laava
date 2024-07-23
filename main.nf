#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { map_reads; make_report } from './modules/local/laava'

NO_FILE = file("$projectDir/bin/NO_FILE")
NO_FILE2 = file("$projectDir/bin/NO_FILE2")

// Unpack the input sample(s) and metadata
def prepareInput(
        seq_reads_file, seq_reads_folder, sample_unique_id, sample_display_name,
        sample_in_metadata
) {
    // Known file extensions
    def SAMPLE_FILE_GLOB = "*.{bam,fastq,fastq.gz,fq,fq.gz}"
    def EXTENSION_REGEX = /\.(bam|fastq|fastq\.gz|fq|fq\.gz)$/

    if (seq_reads_folder) {
        // Multi-sample mode
        def sampleFiles = file("${params.seq_reads_folder}/${SAMPLE_FILE_GLOB}")

        if (params.sample_in_metadata) {
            // Metadata provided - use it to match files
            def csvData = channel
                .fromPath(params.sample_in_metadata)
                .splitCsv(header: true, sep='\t', strip=true)
                .map { row -> [row.sample_unique_id, row.sample_display_name] }
                .toList()
                .val
            def matchedSamples = csvData.collect { sampleId, sampleName ->
                def matchingFiles = sampleFiles.findAll { it.name.contains(sampleId) }
                if (matchingFiles.size() != 1) {
                    error "Error: sample_unique_id '${sampleId}' matches ${matchingFiles.size()} files. Values must be unique."
                }
                [sampleId, sampleName, matchingFiles[0]]
            }
            // Check if all files were matched
            def unmatchedFiles = sampleFiles - matchedSamples.collect { it[2] }
            if (unmatchedFiles) {
                error "Error: The following files were not matched to any sample_unique_id: ${unmatchedFiles.join(', ')}"
            }
            return channel.fromList(matchedSamples)

        } else {
            // No metadata provided - generate sampleId and sampleName from filenames
            return channel.fromList(sampleFiles.collect { sampleFile ->
                def stem = sampleFile.baseName.replaceFirst(EXTENSION_REGEX, '')
                [stem, stem, sampleFile]
            })
        }
    } else if (params.seq_reads_file) {
        // Single-sample mode
        def sampleFile = file(params.seq_reads_file)
        if (!sampleFile.name.matches(/.*${EXTENSION_REGEX}/)) {
            error "Error: The provided sample file '${sampleFile.name}' does not have a supported extension (${SAMPLE_FILE_GLOB})"
        }

        def stem = sampleFile.baseName.replaceFirst(EXTENSION_REGEX, '')
        def sampleId = params.sample_unique_id ?: stem
        def sampleName = params.sample_display_name ?: params.sample_unique_id ?: stem
        return channel.of([sampleId, sampleName, sampleFile])

    } else {
        error "Invalid input parameters. Provide either a sample folder path or a single sample file."
    }
}


workflow laava {
    take:
    seq_reads_file
    seq_reads_folder
    sample_unique_id
    sample_display_name
    sample_in_metadata
    vector_fa
    packaging_fa
    host_fa
    itr_label_1
    itr_label_2
    repcap_name
    helper_name
    lambda_name
    vector_bed
    vector_type
    target_gap_threshold
    max_allowed_outside_vector
    max_allowed_missing_flanking
    flipflop_name
    flipflop_fa

    main:
    // Get a tuple of (ID, name, file) each given sample file and metadata
    sample_channel = prepareInput(
        seq_reads_file, seq_reads_folder, sample_unique_id, sample_display_name,
        sample_in_metadata
    )
    map_reads(
        sample_channel
        .combine( channel.fromPath( vector_fa ) )
        .combine( packaging_fa ? channel.fromPath( packaging_fa ) : channel.of( NO_FILE ) )
        .combine( host_fa ? Channel.fromPath( host_fa ) : Channel.of( NO_FILE2 ) )
        .combine( channel.of( repcap_name ) )
        .combine( channel.of( helper_name ) )
        .combine( channel.of( lambda_name ) )
    )
    make_report(
        map_reads.out.mapped_sam
        .combine( channel.fromPath( vector_bed ) )
        .combine( channel.of( itr_label_1 ) )
        .combine( channel.of( itr_label_2 ) )
        .combine( channel.of( vector_type ) )
        .combine( channel.of( target_gap_threshold ) )
        .combine( channel.of( max_allowed_outside_vector ) )
        .combine( channel.of( max_allowed_missing_flanking ) )
        .combine( channel.of( flipflop_name ) )
        .combine( flipflop_fa ? channel.fromPath( flipflop_fa ) : Channel.of( NO_FILE ) )
    )

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
    laava(
        params.seq_reads_file,
        params.seq_reads_folder,
        params.sample_unique_id,
        params.sample_display_name,
        params.sample_in_metadata,
        params.vector_fa,
        params.packaging_fa,
        params.host_fa,
        params.itr_label_1,
        params.itr_label_2,
        params.repcap_name,
        params.helper_name,
        params.lambda_name,
        params.vector_bed,
        params.vector_type,
        params.target_gap_threshold,
        params.max_allowed_outside_vector,
        params.max_allowed_missing_flanking,
        params.flipflop_name,
        params.flipflop_fa
    )
}
