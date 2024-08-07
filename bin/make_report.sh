#!/bin/bash -ex
sample_name=$1
vector_type=$2
target_gap_threshold=$3
max_allowed_outside_vector=$4
max_allowed_missing_flanking=$5
mapped_reads_sam=$6
annotation_txt=$7
flipflop_name=$8
flipflop_fa=$9

ls -Alh

if [[ $mapped_reads_sam == *.bam ]]; then
    echo "Converting $mapped_reads_sam from BAM to SAM"
    sam_fname="${mapped_reads_sam%%.sam}"
    samtools view -b -o "$sam_fname" "$mapped_reads_sam"
    mapped_reads_sam="$sam_fname"
else
    echo "Reads $mapped_reads_sam appear to be in SAM format"
fi

echo
echo "Starting summarize_AAV_alignment"
summarize_AAV_alignment.py \
    "$mapped_reads_sam" "$annotation_txt" "$sample_name" \
    --target-gap-threshold=$target_gap_threshold \
    --max-allowed-outside-vector=$max_allowed_outside_vector \
    --max-allowed-missing-flanking=$max_allowed_missing_flanking \
    --cpus $(nproc)

echo "Finished summarize_AAV_alignment"
ls -Alh

if [[ -n "$flipflop_name" || -n "$flipflop_fa" ]]; then
    if [ -n "$flipflop_fa" ]; then
        # Use the gives seqs for FF analysis, regardless of FF name
        ff_opt="--flipflop-fasta $flipflop_fa"
    elif [ "$flipflop_name" == "AAV2" ]; then
        # Use the built-in AAV2 sequences
        ff_opt=""
    else
        echo "No built-in sequences for $flipflop_name and none given"
        exit 1
    fi

    echo
    echo "Starting get_flipflop_config"
    get_flipflop_config.py \
        "${sample_name}.tagged.bam" "${sample_name}.per_read.tsv" \
        $ff_opt \
        -o "$sample_name"
    echo "Finished get_flipflop_config"
    ls -Alh
    flipflop_assignments="${sample_name}.flipflop_assignments.tsv"
else
    echo "Skipping flip/flop analysis"
    flipflop_assignments=""
fi

echo
echo "Starting calculate_rdata"
calculate_rdata.R "./${sample_name}" "$annotation_txt" "$sample_name" "$vector_type" "$flipflop_assignments"
echo "Finished calculate_rdata"
ls -Alh

echo
echo "Starting create_report"
create_report.R "./${sample_name}.Rdata"
echo "Finished create_report"
ls -Alh
