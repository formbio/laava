#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { match_metadata_to_files; map_reads; make_report } from './modules/local/laava'

NO_FILE = file("${params.laava_dir}/bin/NO_FILE")
NO_FILE2 = file("${params.laava_dir}/bin/NO_FILE2")

// Unpack the input sample(s) and metadata
def prepare_input(
        seq_reads_file, seq_reads_folder, sample_unique_id, sample_display_name,
        sample_in_metadata
) {
    def SAMPLE_FILE_GLOB = "*.{bam,fastq,fastq.gz,fq,fq.gz}"
    def EXTENSION_REGEX = /\.(bam|fastq|fastq\.gz|fq|fq\.gz)$/

    if (seq_reads_folder) {
        // Multi-sample mode
        if (sample_in_metadata) {
            // TSV provided - load it in a separate process
            return match_metadata_to_files(file(sample_in_metadata), file(seq_reads_folder))
                .splitCsv(sep: '\t')
                .map { row -> [row[0], row[1], file("${seq_reads_folder}/" + row[2])] }

        } else {
            // No TSV provided - generate sample_id and sample_name from filenames
            //def found_files = file("${sample_folder}/*.{bam,fastq,fastq.gz,fq,fq.gz}")
            def found_files = file("${seq_reads_folder}/${SAMPLE_FILE_GLOB}")
            return channel.fromList(found_files.collect { seqfile ->
                def stem = seqfile.baseName.replaceFirst(EXTENSION_REGEX, '')
                [stem, stem, seqfile]
            })

        }
    } else if (seq_reads_file) {
        // Single-sample mode
        def seq_file = file(seq_reads_file)
        if (!seq_file.exists()) {
            error "Error: The provided sample file '${seq_reads_file}' does not exist."
        }
        if (!seq_file.name.matches(/.*${EXTENSION_REGEX}/)) {
            error "Error: The provided sample file '${seq_file.name}' does not have a supported extension (bam, fastq, fastq.gz, fq, fq.gz)"
        }

        def stem = seq_file.baseName.replaceFirst(EXTENSION_REGEX, '')
        def sample_id = sample_unique_id ?: stem
        def sample_name = sample_display_name ?: sample_unique_id ?: stem
        return channel.of([sample_id, sample_name, seq_file])

    } else {
        error "Invalid input parameters. Provide either a sample folder, a TSV file with sample folder, or a single sample file."
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
    mitr_label
    repcap_name
    helper_name
    lambda_name
    vector_bed
    vector_type
    target_gap_threshold
    max_allowed_outside_vector
    max_allowed_missing_flanking
    min_supp_joint_coverage
    flipflop_name
    flipflop_fa

    main:
    // Get a tuple of (ID, name, file) each given sample file and metadata
    sample_channel = prepare_input(
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
        .combine( channel.of( mitr_label ) )
        .combine( channel.of( vector_type ) )
        .combine( channel.of( target_gap_threshold ) )
        .combine( channel.of( max_allowed_outside_vector ) )
        .combine( channel.of( max_allowed_missing_flanking ) )
        .combine( channel.of( min_supp_joint_coverage ) )
        .combine( channel.of( flipflop_name ) )
        .combine( flipflop_fa ? channel.fromPath( flipflop_fa ) : Channel.of( NO_FILE ) )
    )

    emit:
    mapped_sam = map_reads.out.mapped_sam
    mapped_bam = map_reads.out.mapped_bam
    metadata_out_tsv = make_report.out.metadata_tsv
    alignments_tsv = make_report.out.alignments_tsv
    per_read_tsv = make_report.out.per_read_tsv
    nonmatch_tsv = make_report.out.nonmatch_tsv
    agg_ref_type_tsv = make_report.out.agg_ref_type_tsv
    agg_subtype_tsv = make_report.out.agg_subtype_tsv
    agg_flipflop_tsv = make_report.out.agg_flipflop_tsv
    tagged_bam = make_report.out.tagged_bam
    subtype_bams = make_report.out.subtype_bams
    subtype_bais = make_report.out.subtype_bais
    flipflop_bams = make_report.out.flipflop_bams
    flipflop_tsv = make_report.out.flipflop_tsv
    sample_id = make_report.out.sample_id
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
        params.mitr_label,
        params.repcap_name,
        params.helper_name,
        params.lambda_name,
        params.vector_bed,
        params.vector_type,
        params.target_gap_threshold,
        params.max_allowed_outside_vector,
        params.max_allowed_missing_flanking,
        params.min_supp_joint_coverage,
        params.flipflop_name,
        params.flipflop_fa
    )
}
