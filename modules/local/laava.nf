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
          path("${sample_id}.reference_names.tsv"),
          path("${sample_id}.sort_by_name.sam"), emit: mapped_sam
    tuple val(sample_id),
          val(sample_name),
          path("${sample_id}.sort_by_pos.bam"),
          path("${sample_id}.sort_by_pos.bam.bai"), emit: mapped_bam

    script:
    // Hack for optional inputs
    def packaging_fa_opt = packaging_fa.name != "NO_FILE" ? "--packaging \"$packaging_fa\"" : ""
    def host_fa_opt = host_fa.name != "NO_FILE2" ? "--host \"$host_fa\"" : ""
    def repcap_name_opt = repcap_name ? "--repcap-name \"$repcap_name\"" : ""
    def helper_name_opt = helper_name ? "--helper-name \"$helper_name\"" : ""
    def lambda_name_opt = lambda_name ? "--lambda-name \"$lambda_name\"" : ""
    def packaging_fa_path = packaging_fa.name != "NO_FILE" ? "$packaging_fa" : ""
    def host_fa_path = host_fa.name != "NO_FILE2" ? "$host_fa" : ""
    """
    get_reference_names.py "${vector_fa}" ${packaging_fa_opt} ${host_fa_opt} \\
        ${repcap_name_opt} ${helper_name_opt} ${lambda_name_opt} \\
        -o "${sample_id}.reference_names.tsv"
    map_reads.sh ${sample_id} "${reads}" "${vector_fa}" \\
        "${packaging_fa_path}" "${host_fa_path}"
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
    // summary tables
    path("${sample_id}.metadata.tsv"), emit: metadata_tsv
    path("${sample_id}.alignments.tsv.gz"), emit: alignments_tsv
    path("${sample_id}.per_read.tsv.gz"), emit: per_read_tsv
    path("${sample_id}.nonmatch.tsv.gz"), emit: nonmatch_tsv
    path("${sample_id}.flipflop.tsv.gz"), emit: flipflop_tsv, optional: true
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
    write_sample_metadata.py "${sample_id}" "${sample_name}" "${mapped_reads}" \\
        -o "${sample_id}.metadata.tsv"
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
