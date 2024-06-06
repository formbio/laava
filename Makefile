# Build Docker images
all: laava laava_dev

.PHONY: clean laava laava_dev diffcheck
clean:
	rm -fv .nextflow.log*
	rm -fv test/build/*
	rm -rf workflow-outputs/*

laava laava_dev: %: %.dockerfile laava.conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .

diffcheck:
	diff test/build-snapshot/sc.per_read.csv workflow-outputs/output/sc.per_read.csv && echo "OK"
	diff test/build-snapshot/ss.per_read.csv workflow-outputs/output/ss.per_read.csv && echo "OK"
	diff test/build-snapshot/ss.flipflop.tsv workflow-outputs/output/ss.flipflop.tsv && echo "OK"
