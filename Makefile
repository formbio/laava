# Build Docker images and test Nextflow

# Nextflow workflow output directory
wf_out_dir := workflow-outputs/output

all: laava laava_dev

.PHONY: clean laava laava_dev sc ss diffcheck-sc diffcheck-ss
clean:
	rm -fv .nextflow.log*
	rm -fv test/build/*
	rm -rf workflow-outputs/*

laava laava_dev: %: %.dockerfile laava.conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .


sc: params-local-sc-no-ff.json
	nextflow run -profile local main.nf -params-file $<

ss: params-local-ss-with-ff.json
	nextflow run -profile local main.nf -params-file $<

min: params-local-no-file-sc.json
	nextflow run -profile local main.nf -params-file $<

diffcheck-sc: $(wf_out_dir)/sc.subsample005.bam.per_read.tsv
	diff test/build-snapshot/sc.per_read.tsv $< && echo "OK"

diffcheck-ss:  $(wf_out_dir)/ss.subsample005.bam.per_read.tsv $(wf_out_dir)/ss.subsample005.bam.flipflop.tsv
	diff test/build-snapshot/ss.per_read.tsv $< && echo "OK"
	diff test/build-snapshot/ss.flipflop.tsv $(lastword $^) && echo "OK"
