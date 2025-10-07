#!/bin/bash -ex
# Implement the 'make_report' Nextflow process.
sample_id=$1
sample_name=$2
version=$3
reference_names=$4
mapped_reads=$5
vector_annotation=$6
itr_label_1=$7
itr_label_2=$8
mitr_label=$9
vector_type=${10}
target_gap_threshold=${11}
max_allowed_outside_vector=${12}
max_allowed_missing_flanking=${13}
min_supp_joint_coverage=${14}
flipflop_name=${15}
flipflop_fa=${16}
out_dir=${17}


ls -Alh

python "$(dirname $0)/write_sample_metadata.py" "$sample_id" "$sample_name" "$mapped_reads" \
    -v "$version" -o "$out_dir/${sample_id}.metadata.tsv"

ls -Alh

if [ "$vector_type" == "unspecified" ]; then
    vector_type=$(python "$(dirname $0)/guess_vector_type_length.py" "$vector_annotation" \
                  "$itr_label_1" "$itr_label_2" "$mitr_label")
    echo "Inferred vector_type: $vector_type"
fi

echo
echo "Starting summarize_alignment"
python "$(dirname $0)/summarize_alignment.py" \
    "$mapped_reads" "$vector_annotation" "${reference_names}" \
    "${itr_label_1}" "${itr_label_2}" "${mitr_label}" \
    --output-prefix="$out_dir/$sample_id" \
    --sample-id="$sample_id" \
    --vector-type="$vector_type" \
    --target-gap-threshold="$target_gap_threshold" \
    --max-allowed-outside-vector="$max_allowed_outside_vector" \
    --max-allowed-missing-flanking="$max_allowed_missing_flanking" \
    --min-supp-joint-coverage="$min_supp_joint_coverage" \
    --cpus $(if command -v nproc >/dev/null 2>&1; then nproc; else sysctl -n hw.ncpu 2>/dev/null || echo 1; fi)

echo "Finished summarize_alignment"
ls -Alh
: <<'TEST'   
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
    
    if [ "${vector_type}" == "sc" ]; then
       m_itr_start=$(cat "$vector_annotation" | grep "${mitr_label}" | cut -f2)
       wt_start=$(cat "$vector_annotation" | grep "${itr_label_1}" | cut -f2)
       if [ "$wt_start" -gt "$m_itr_start" ]; then
          orientation="left"
       elif [ "$m_itr_start" -gt "$wt_start" ]; then
          orientation="right"
       fi
    elif [ "${vector_type}" == "ss" ]; then
       orientation="left"
    fi
 
    python "$(dirname $0)/get_flipflop_config.py" \
        "$out_dir/${sample_id}.tagged.bam" "$out_dir/${sample_id}.per_read.tsv.gz" \
        "$vector_type" \
        "$orientation" \
        $ff_opt \
        -o "$out_dir/$sample_id"
    echo "Finished get_flipflop_config"
    ls -Alh
else
    echo "Skipping flip/flop analysis"
fi

echo "Starting aggregate_tables"
python "$(dirname $0)/aggregate_tables.py" --path-prefix "$out_dir/$sample_id"
echo "Finished aggregate_tables"

ls -Alh

echo
echo "Starting create_report"
if [ -n "$CONDA_PREFIX" ]; then
    # Use conda's R with proper library path
    export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
    $CONDA_PREFIX/bin/Rscript "$(dirname $0)/create_report.R" "$out_dir/$sample_id" "$vector_type" \
        $(python "$(dirname $0)/emit_target_coords.py" "${vector_annotation}" \
            "${itr_label_1}" "${itr_label_2}" "${mitr_label}")
else
    # Fallback to system R
    create_report.R "$out_dir/$sample_id" "$vector_type" \
        $(python "$(dirname $0)/emit_target_coords.py" "${vector_annotation}" \
            "${itr_label_1}" "${itr_label_2}" "${mitr_label}")
fi
echo "Finished create_report"

ls -Alh
TEST