#!/bin/bash -e
# Extracted from bin/map_reads.sh
vector_fa=$1
packaging_fa=$2
host_fa=$3
repcap_name=$4
out_fname=$5


function get_reference_names {
    fasta=$1
    typename=$2
    grep '^>' "$fasta" \
        | cut -f1 -d\  | cut -f2 -d\> \
        | awk -v var="$typename" '{
            print $1,var
        }'
}

function relabel_repcap {
    tsv=$1
    rc_name=$2
    cat "$tsv" | awk -v var="$rc_name" '{
        if ($1==var) {
            print $1,"repcap"
        } else {
            print $1,$2
        }
    }'
}


# Consolidate inputs, skipping optional files
get_reference_names "$vector_fa" "vector" > reference_names.tsv
if [ -e "$packaging_fa" ]; then
    get_reference_names "$packaging_fa" "helper" > _tmp.tsv
    # Re-label 'repcap' within packaging_fa
    if [ -n "$repcap_name" ]; then
        relabel_repcap _tmp.tsv "$repcap_name" \
            >> reference_names.tsv
    else
        cat _tmp.tsv >> reference_names.tsv
    fi
fi
if [ -e "$host_fa" ]; then
    get_reference_names "$host_fa" "host" >> reference_names.tsv
fi

mv reference_names.tsv $out_fname
