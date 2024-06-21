#!/bin/bash -ex
sample_name=$1
reads=$2
vector_fa=$3
packaging_fa=$4
host_fa=$5
repcap_name=$6

# NOTE: the sequence IDs should be free of blank spaces and symbols. Stick with numbers,
# alphabet letters, and _ and -. If necessary, rename the sequence IDs in the combined
# fasta file.

# STEPS:
# - cat the reference fasta sequences
# - [enh] fix reference fasta sequence names
# - align reads to references OUTPUT
# - BAM, aligned to all references and sorted by read name (?!) SOMEWHERE
# - extract all references' names and lengths (.genome, via faidx)
# DOWNSTREAM:
# - take vector annotations (BED or TXT)
# - append reference names and lengths to vector annotations TXT

ls -Alh

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
cp "$vector_fa" all_refs.fa
get_reference_names "$vector_fa" "vector" > reference_names.tsv
if [ -e "$packaging_fa" ]; then
    cat "$packaging_fa" >> all_refs.fa
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
    cat "$host_fa" >> all_refs.fa
    get_reference_names "$host_fa" "host" >> reference_names.tsv
fi

# Logging
cat reference_names.tsv
echo
grep '^>' all_refs.fa
echo

if [[ $reads == *.bam ]]; then
    echo "Converting $reads from BAM to FASTQ"
    samtools fastq -n -0 reads.fq "$reads"
    reads_fn=reads.fq
else
    # NB: minimap2 handles .gz automatically
    reads_fn="$reads"
fi

threads=$(nproc)
minimap2 --eqx -a --secondary=no -t $threads all_refs.fa "$reads_fn" > tmp.mapped.sam
# Sort the mapped reads by name
samtools sort -@ $threads -n -O SAM -o "$sample_name.sort_by_name.sam" tmp.mapped.sam

# Make a position-sorted BAM output file for other downstream consumers
out_bam="$sample_name.sort_by_pos.bam"
# Drop unmapped reads
samtools view -@ $threads --fast -F 4 -o tmp.sorted.bam tmp.mapped.sam
samtools sort -@ $threads -o "$out_bam" tmp.sorted.bam
samtools index "$out_bam"

ls -Alh
