#!/bin/bash -ex
# Implement the 'make_report' Nextflow process.
sample_id=$1
sample_name=$2
reference_names=$3
mapped_reads=$4
vector_annotation=$5
itr_label_1=$6
itr_label_2=$7
mitr_label=$8
vector_type=$9
target_gap_threshold=${10}
max_allowed_outside_vector=${11}
max_allowed_missing_flanking=${12}
flipflop_name=${13}
flipflop_fa=${14}

# Special filename for the intermediate annotation format
annotation_txt="annotation.txt"

ls -Alh

write_sample_metadata.py "${sample_id}" "${sample_name}" "${mapped_reads}" \
    -o "${sample_id}.metadata.tsv"

ls -Alh

prepare_annotation.py "${vector_annotation}" "${reference_names}" \
    "${itr_label_1}" "${itr_label_2}" "${mitr_label}" \
    -o "$annotation_txt"


ls -Alh

if [[ $mapped_reads == *.bam ]]; then
    echo "Converting $mapped_reads from BAM to SAM"
    sam_fname="${mapped_reads%%.sam}"
    samtools view -b -o "$sam_fname" "$mapped_reads"
    mapped_reads="$sam_fname"
else
    echo "Reads $mapped_reads appear to be in SAM format"
fi

if [ "$vector_type" == "unspecified" ]; then
    vector_type=$(guess_vector_type_length.py "$vector_annotation" \
                  "${itr_label_1}" "${itr_label_2}" "${mitr_label}")
    echo "Inferred vector_type: $vector_type"
fi

echo
echo "Starting summarize_alignment"
summarize_alignment.py \
    "$mapped_reads" "$vector_annotation" "${reference_names}" \
    "${itr_label_1}" "${itr_label_2}" "${mitr_label}" \
    --output-prefix="$sample_id" \
    --sample-id="$sample_id" \
    --vector-type="$vector_type" \
    --target-gap-threshold="$target_gap_threshold" \
    --max-allowed-outside-vector="$max_allowed_outside_vector" \
    --max-allowed-missing-flanking="$max_allowed_missing_flanking" \
    --cpus $(nproc)

echo "Finished summarize_alignment"
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
        "${sample_id}.tagged.bam" "${sample_id}.per_read.tsv.gz" \
        $ff_opt \
        -o "$sample_id"
    echo "Finished get_flipflop_config"
    ls -Alh
else
    echo "Skipping flip/flop analysis"
fi

echo "Starting aggregate_tables"
aggregate_tables.py --path-prefix "${sample_id}"
echo "Finished aggregate_tables"

ls -Alh

echo
echo "Starting create_report"
create_report.R "./${sample_id}" "$vector_type" "$annotation_txt"
echo "Finished create_report"

ls -Alh
