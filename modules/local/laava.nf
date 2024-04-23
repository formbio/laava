process map_reads() {
    publishDir "$params.output", mode: "copy"

    input:
    val sample_name
    path reads
    path vector_fa
    path helper_fa
    path repcap_fa
    path host_fa

    output:
    path("${sample_name}.sort_by_name.sam"), emit: mapped_reads
    path("reference_names.tsv"), emit: reference_names
    val "${sample_name}", emit: sample_name

    script:
    // Hack for optional inputs
    def helper_fa_path = helper_fa.name != "NO_FILE" ? "$helper_fa" : ""
    def repcap_fa_path = repcap_fa.name != "NO_FILE" ? "$repcap_fa" : ""
    def host_fa_path = host_fa.name != "NO_FILE" ? "$host_fa" : ""
    """
    map_reads.sh ${sample_name} "${reads}" "${vector_fa}" \\
        "${helper_fa_path}" "${repcap_fa_path}" "${host_fa_path}"
    """
}


process make_report() {
    publishDir "$params.output", mode: "copy"

    input:
    path vector_annotation
    path reference_names
    val sample_name
    val target_gap_threshold
    val max_allowed_outside_vector
    val max_allowed_missing_flanking
    path mapped_reads
    val flipflop_name
    path flipflop_fa

    output:
    // summarize alignment
    path("${sample_name}.per_read.csv"), emit: per_read_csv
    path("${sample_name}.summary.csv"), emit: summary_csv
    path("${sample_name}.nonmatch_stat.csv.gz"), emit: nonmatch_stat_csvgz
    path("${sample_name}.tagged.bam"), emit: tagged_bam
    path("${sample_name}.*.tagged.sorted.bam"), emit: subtype_bams
    path("${sample_name}.*.tagged.sorted.bam.bai"), emit: subtype_bais
    // flip-flop
    path("${sample_name}.flipflop_assignments.txt"), emit: flipflop_assignments_txt, optional: true
    path("${sample_name}.*-flipflop.bam"), emit: flipflop_bams
    // report
    path("${sample_name}.alignments.tsv"), emit: alignments_tsv
    path("${sample_name}.readsummary.tsv"), emit: readsummary_tsv
    path("${sample_name}.sequence-error.tsv"), emit: sequence_error_tsv
    path("${sample_name}.flipflop.tsv"), emit: flipflop_tsv, optional: true
    path("${sample_name}.Rdata"), emit: rdata, optional: true
    path("${sample_name}_AAV_report.html"), emit: aav_report_html
    path("${sample_name}_AAV_report.pdf"), emit: aav_report_pdf

    script:
    def ff_fa_path = flipflop_fa.name != "NO_FILE" ? "$flipflop_fa" : ""
    """
    prepare_annotation.py "${vector_annotation}" "${reference_names}" -o annotation.txt
    make_report.sh "${sample_name}" \\
        $target_gap_threshold $max_allowed_outside_vector $max_allowed_missing_flanking \\
        "${mapped_reads}" annotation.txt \\
        "${flipflop_name}" "${ff_fa_path}"
    """
}
