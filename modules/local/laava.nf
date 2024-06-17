process map_reads() {
   label 'laava'
    publishDir "$params.output", mode: "copy"

    input:
    tuple val(sample_name),
          path(reads),
          path(vector_fa),
          path(packaging_fa),
          path(host_fa),
          val(repcap_name)

    output:
    tuple val(sample_name),
          path("reference_names.tsv"),
          path("${sample_name}.sort_by_name.sam"), emit: mapped_reads
    tuple val(sample_name),
          path("${sample_name}.bam"), emit: bam
    script:
    // Hack for optional inputs
    def packaging_fa_path = packaging_fa.name != "NO_FILE" ? "$packaging_fa" : ""
    def host_fa_path = host_fa.name != "NO_FILE" ? "$host_fa" : ""
    """
    map_reads.sh ${sample_name} "${reads}" "${vector_fa}" \\
        "${packaging_fa_path}" "${host_fa_path}" "${repcap_name}" 
    """
}
process hostgenect {
  label 'laavasupp'
  publishDir "${params.output}", mode: 'copy'
  input:
  tuple val(sid),path(sbam),path(exonbed)
  output:
  path("${sid}.bedtools.cov.txt")
  script:
  """
  hostgenect.sh -p ${sid} -e ${exonbed} ${sbam}
  """
}
process bamqc {
  label 'laavasupp'
  publishDir "$params.output/qc", mode: 'copy'
  input:
  tuple val(sid),path(bam)
  output:
  path("${sid}*"), emit: qc optional true
  script:
  """
  bamqc.sh -b ${bam} -p ${sid} 
  """
}
process make_report() {
   label 'laava'
   publishDir "$params.output", mode: "copy"

    input:
    tuple val(sample_name),
          path(reference_names),
          path(mapped_reads),
          path(vector_annotation),
          val(vector_type),
          val(target_gap_threshold),
          val(max_allowed_outside_vector),
          val(max_allowed_missing_flanking),
          val(flipflop_name),
          path(flipflop_fa)

    output:
    // summarize alignment
    path("${sample_name}.per_read.tsv"), emit: per_read_tsv
    path("${sample_name}.summary.tsv"), emit: summary_tsv
    path("${sample_name}.nonmatch_stat.tsv.gz"), emit: nonmatch_stat_tsvgz
    path("${sample_name}.tagged.bam"), emit: tagged_bam
    path("${sample_name}.*.tagged.sorted.bam"), emit: subtype_bams
    path("${sample_name}.*.tagged.sorted.bam.bai"), emit: subtype_bais
    // flip-flop
    path("${sample_name}.flipflop_assignments.tsv"), emit: flipflop_assignments_tsv, optional: true
    path("${sample_name}.*-flipflop.bam"), emit: flipflop_bams, optional: true
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
    make_report.sh \\
        "${sample_name}" \\
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
