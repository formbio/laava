# Build Docker images and test Nextflow

# Nextflow workflow output directory
wf_out_dir := workflow-outputs/output
snapshot_dir := test/build-snapshot

all: laava laava_dev

.PHONY: clean laava laava_dev sc ss diffcheck-sc diffcheck-ss
clean:
	rm -f .nextflow.log*
	rm -fr .nextflow/*
	rm -fr workflow-outputs/*


laava laava_dev: %: %.dockerfile laava.conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .


sc: params-local-sc-no-ff.json
	nextflow run -profile local main.nf -params-file $<

ss: params-local-ss-with-ff.json
	nextflow run -profile local main.nf -params-file $<

min: params-local-no-file-sc.json
	nextflow run -profile local main.nf -params-file $<

folder: params-local-sc-folder.json
	nextflow run -profile local main.nf -params-file $<


diffcheck-sc: $(wf_out_dir)/sc.subsample005.per_read.tsv
	diff $(snapshot_dir)/sc.per_read.tsv $< && echo "OK"

diffcheck-ss: $(wf_out_dir)/ss.subsample005.per_read.tsv $(wf_out_dir)/ss.subsample005.flipflop.tsv
	#diff $(snapshot_dir)/ss.per_read.tsv $< && echo "OK"
	diff $(snapshot_dir)/ss.flipflop.tsv $(lastword $^) && echo "OK"
