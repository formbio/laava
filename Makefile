# Build Docker images and test Nextflow

# Nextflow workflow output directory
wf_out_dir := workflow-outputs/output

# Local test directories
in_dir = test/samples
out_dir = test/build
fa_dir = test/fasta

# Form Bio workflow deployment
formbio_org := form-bio-solutions
formbio_project := aav-qc-workshop
# Avoid uploading the local test dir; it's > the upload size limit
tmp_stash_dir := /tmp/laava-deploy-test
docker_repo := ghcr.io/formbio

all: laava laava_dev sc ss min folder test lint

.PHONY: clean laava laava_dev formbio sc ss min folder test test-local lint lint-from-docker
clean:
	rm -f .nextflow.log*
	rm -fr .nextflow/*
	rm -fr workflow-outputs/*
	rm -frv $(out_dir)/*

# Build & deploy

laava: Dockerfile conda_env.yml
	docker build -t $(docker_repo)/$@:dev -f $< .

laava_dev: laava_dev.dockerfile conda_env.yml
	docker build -t $(docker_repo)/$@:latest -f $< .

formbio: clean
	mv test/ "$(tmp_stash_dir)"
	formbio workflow upload \
		--org "$(formbio_org)" --project "$(formbio_project)" \
		--env prod --visibility PROJECT \
		--version dev --repo . --workflow laava
	mv "$(tmp_stash_dir)" test


# Test Nextflow pipeline execution

sc ss min folder: %: params-local-%.json
	nextflow run -profile local main.nf -params-file $<

# Test local execution and outputs

test: lint-from-docker laava_dev
	docker run --rm -v $(CURDIR):/data -w /data -it $(docker_repo)/laava_dev:latest \
		make -B -C test test

lint-from-docker: laava_dev
	docker run --rm -v $(CURDIR):/data -w /data -it $(docker_repo)/laava_dev:latest \
		make lint

# Test local execution without Docker (using conda)

test-local: lint
	cd test && make test-local

lint: lint-r lint-python

lint-python:
	ruff check --isolated --output-format=github --no-cache src/

lint-r:
	 Rscript -e 'library(lintr); options(lintr.error_on_lint=TRUE); lint_dir(".", linters=linters_with_tags("correctness"))'
