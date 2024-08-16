#!/bin/bash -e
snapshot_dir=$1
build_dir=$2
shift 2
tsv_paths=$@

for fpath in $tsv_paths; do
    fname=$(basename $fpath)
    expect_path="${snapshot_dir}/${fname}"
    observ_path="${build_dir}/${fname}"
    test $(cat "$expect_path" | wc -l) -eq $(cat "$observ_path" | wc -l)
	diff "$expect_path" "$observ_path"
done
