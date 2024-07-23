process map_reads() {
    publishDir "$params.output", mode: "copy"

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
          path("reference_names.tsv"),
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
    map_reads.sh ${sample_id} "${reads}" "${vector_fa}" \\
        "${packaging_fa_path}" "${host_fa_path}" \\
        "${repcap_name}" "${helper_name}" "${lambda_name}"
    """
}


process make_report() {
    publishDir "$params.output", mode: "copy"

    input:
    tuple val(sample_id),
          val(sample_name),
          path(reference_names),
          path(mapped_reads),
          path(vector_annotation),
          val(itr_label_1),
          val(itr_label_2),
          val(vector_type),
          val(target_gap_threshold),
          val(max_allowed_outside_vector),
          val(max_allowed_missing_flanking),
          val(flipflop_name),
          path(flipflop_fa)

    output:
    // summarize alignment
    path("${sample_id}.per_read.tsv"), emit: per_read_tsv
    path("${sample_id}.summary.tsv"), emit: summary_tsv
    path("${sample_id}.nonmatch_stat.tsv.gz"), emit: nonmatch_stat_tsvgz
    path("${sample_id}.tagged.bam"), emit: tagged_bam
    path("${sample_id}.*.tagged.sorted.bam"), emit: subtype_bams
    path("${sample_id}.*.tagged.sorted.bam.bai"), emit: subtype_bais
    // flip-flop
    path("${sample_id}.flipflop_assignments.tsv"), emit: flipflop_assignments_tsv, optional: true
    path("${sample_id}.*-flipflop.bam"), emit: flipflop_bams, optional: true
    // intermediate data
    path("${sample_id}.alignments.tsv"), emit: alignments_tsv
    path("${sample_id}.readsummary.tsv"), emit: readsummary_tsv
    path("${sample_id}.sequence-error.tsv"), emit: sequence_error_tsv
    path("${sample_id}.flipflop.tsv"), emit: flipflop_tsv, optional: true
    path("${sample_id}.Rdata"), emit: rdata, optional: true
    // report
    path("${sample_id}_AAV_report.html"), emit: aav_report_html
    path("${sample_id}_AAV_report.pdf"), emit: aav_report_pdf

    script:
    def ff_fa_path = flipflop_fa.name != "NO_FILE" ? "$flipflop_fa" : ""
    """
    prepare_annotation.py "${vector_annotation}" "${reference_names}" \\
        "${itr_label_1}" "${itr_label_2}" \\
        -o annotation.txt
    make_report.sh \\
        "${sample_id}" \\
        $vector_type \\
        $target_gap_threshold \\
        $max_allowed_outside_vector \\
        $max_allowed_missing_flanking \\
        "${mapped_reads}" \\
        annotation.txt \\
        "${flipflop_name}" \\
        "${ff_fa_path}"
    """
}
