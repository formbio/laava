#!/bin/bash -ex
sample_id=$1
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

if [ "$vector_type" == "unspecified" ]; then
    #vector_name=$(samtools view -H ss.tagged.bam | head -n 2 | tail -n1 | cut -f2 | cut -d: -f2)
    vector_name=$(head -n1 "$annotation_txt" | cut -d';' -f 1 | cut -d'=' -f2)
    vector_type=$(guess_vector_type.py "$mapped_reads_sam" -v "$vector_name")
    echo "Inferred vector_type: $vector_type"
fi

echo
echo "Starting summarize_alignment"
summarize_alignment.py \
    "$mapped_reads_sam" "$annotation_txt" "$sample_id" \
    --sample-id="$sample_id" \
    --vector-type="$vector_type" \
    --target-gap-threshold=$target_gap_threshold \
    --max-allowed-outside-vector=$max_allowed_outside_vector \
    --max-allowed-missing-flanking=$max_allowed_missing_flanking \
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
    flipflop_tsv="${sample_id}.flipflop.tsv.gz"
else
    echo "Skipping flip/flop analysis"
    flipflop_tsv=""
fi

echo
echo "Starting create_report"
create_report.R "./${sample_id}" "$sample_id" "$vector_type" "$annotation_txt" "$flipflop_tsv"
echo "Finished create_report"
ls -Alh
