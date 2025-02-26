process match_metadata_to_files {
    input:
    path sample_in_metadata
    path sample_folder

    output:
    path("metadata_with_paths.tsv")

    script:
    """
    match_metadata_to_files.py ${sample_in_metadata} ${sample_folder} \\
        > metadata_with_paths.tsv
    """
}


process map_reads() {
    container "${params.container_repo}/laava${params.container_version}"
    publishDir "$params.output/alignment", mode: "copy"

    input:
    tuple val(sample_id),
          val(sample_name),
          path(reads),
          path(vector_fa),
          path(packaging_fa),
          path(host_fa),
          val(repcap_name),
          val(helper_name),
          val(lambda_name)

    output:
    tuple val(sample_id),
          val(sample_name),
          path("${sample_id}.reference_names.tsv"),
          path("${sample_id}.sort_by_name.sam"), emit: mapped_sam
    tuple val(sample_id),
          val(sample_name),
          path("${sample_id}.sort_by_pos.bam"),
          path("${sample_id}.sort_by_pos.bam.bai"), emit: mapped_bam

    script:
    // Hack for optional inputs
    def packaging_fa_path = packaging_fa.name != "NO_FILE" ? "$packaging_fa" : ""
    def host_fa_path = host_fa.name != "NO_FILE2" ? "$host_fa" : ""
    """
    map_reads.sh ${sample_id} "${reads}" \\
        "${vector_fa}" "${packaging_fa}" "${host_fa}" \\
        "${repcap_name}" "${helper_name}" "${lambda_name}"
    """
}


process make_report() {
    container "${params.container_repo}/laava${params.container_version}"
    publishDir "$params.output/report", mode: "copy"

    input:
    tuple val(sample_id),
          val(sample_name),
          path(reference_names),
          path(mapped_reads),
          path(vector_annotation),
          val(itr_label_1),
          val(itr_label_2),
          val(mitr_label),
          val(vector_type),
          val(target_gap_threshold),
          val(max_allowed_outside_vector),
          val(max_allowed_missing_flanking),
          val(min_supp_joint_coverage),
          val(flipflop_name),
          path(flipflop_fa)

    output:
    val(sample_id), emit: sample_id
    // summary tables
    path("${sample_id}.metadata.tsv"), emit: metadata_tsv
    path("${sample_id}.alignments.tsv.gz"), emit: alignments_tsv
    path("${sample_id}.per_read.tsv.gz"), emit: per_read_tsv
    path("${sample_id}.nonmatch.tsv.gz"), emit: nonmatch_tsv
    path("${sample_id}.agg_ref_type.tsv"), emit: agg_ref_type_tsv
    path("${sample_id}.agg_subtype.tsv"), emit: agg_subtype_tsv
    path("${sample_id}.flipflop.tsv.gz"), emit: flipflop_tsv, optional: true
    path("${sample_id}.agg_flipflop.tsv"), emit: agg_flipflop_tsv, optional:true

    // intermediate data
    path("${sample_id}.tagged.bam"), emit: tagged_bam
    path("${sample_id}.*.tagged.sorted.bam"), emit: subtype_bams
    path("${sample_id}.*.tagged.sorted.bam.bai"), emit: subtype_bais
    path("${sample_id}.flipflop-*.bam"), emit: flipflop_bams, optional: true
    // report
    path("${sample_id}_AAV_report.html"), emit: aav_report_html
    path("${sample_id}_AAV_report.pdf"), emit: aav_report_pdf

    script:
    def ff_fa_path = flipflop_fa.name != "NO_FILE" ? "$flipflop_fa" : ""
    """
    make_report.sh \\
        "${sample_id}" \\
        "${sample_name}" \\
        "${workflow.manifest.version}" \\
        "${reference_names}" \\
        "${mapped_reads}" \\
        "${vector_annotation}" \\
        "${itr_label_1}" \\
        "${itr_label_2}" \\
        "${mitr_label}" \\
        "${vector_type}" \\
        "${target_gap_threshold}" \\
        "${max_allowed_outside_vector}" \\
        "${max_allowed_missing_flanking}" \\
        "${min_supp_joint_coverage}" \\
        "${flipflop_name}" \\
        "${ff_fa_path}"
    """
}
