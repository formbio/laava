# Build Docker images and test Nextflow

# Nextflow workflow output directory
wf_out_dir := workflow-outputs/output
snapshot_dir := test/build-snapshot

# Form Bio workflow deployment
formbio_org := form-bio-solutions
formbio_project := aav-qc-workshop
# Avoid uploading the local test dir; it's > the upload size limit
tmp_stash_dir := /tmp/laava-deploy-test

all: laava laava_dev sc ss min folder

.PHONY: clean laava laava_dev sc ss min folder diffcheck-sc diffcheck-ss formbio
clean:
	rm -f .nextflow.log*
	rm -fr .nextflow/*
	rm -fr workflow-outputs/*


laava: Dockerfile conda_env.yml
	docker build -t ghcr.io/formbio/$@:dev -f $< .

laava_dev: laava_dev.dockerfile conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .


sc ss min folder: %: params-local-%.json
	nextflow run -profile local main.nf -params-file $<


diffcheck-sc: $(wf_out_dir)/sc.subsample005.per_read.tsv
	diff $(snapshot_dir)/sc.per_read.tsv $< && echo "OK"

diffcheck-ss: $(wf_out_dir)/ss.subsample005.per_read.tsv $(wf_out_dir)/ss.subsample005.flipflop.tsv
	#diff $(snapshot_dir)/ss.per_read.tsv $< && echo "OK"
	diff $(snapshot_dir)/ss.flipflop.tsv $(lastword $^) && echo "OK"


formbio: clean
	mv test/ "$(tmp_stash_dir)"
	formbio workflow upload \
		--org "$(formbio_org)" --project "$(formbio_project)" \
		--env prod --visibility PROJECT \
		--version dev --repo . --workflow laava
	mv "$(tmp_stash_dir)" test
