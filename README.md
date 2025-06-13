[![Test Status](https://github.com/formbio/laava/actions/workflows/tests.yaml/badge.svg)](https://github.com/formbio/laava/actions/workflows/tests.yaml)

# LAAVA: Long-read AAV Analysis

A software package and executable bioinformatics workflow for the analysis of recombinant adeno-associated virus (rAAV) products by PacBio long-read sequencing.

* For a full explanation of the methods and example results on public PacBio datasets,
  see the preprint paper on bioRxiv:
  [Standardized Nomenclature and Reporting for PacBio HiFi Sequencing and Analysis of rAAV Gene Therapy Vectors](https://www.biorxiv.org/content/10.1101/2024.05.07.592296v1)
* For a summary of technical methods, AAV type/subtype definitions, and interpretation,
  see: [Design and definitions](https://github.com/formbio/laava/wiki/Design-and-definitions)
* For answers to frequently asked questions, see [the FAQ](https://github.com/formbio/laava/wiki/Frequently-Asked-Questions-(FAQ))

## Installation & Usage

LAAVA can be used as an end-to-end Nextflow workflow, an interactive Docker container,
or individual scripts in this codebase.

### Standard Nextflow (recommended)

This code can be run as a standard [Nextflow](https://www.nextflow.io/) workflow.
When run this way, the workflow will automatically pull in the analysis scripts and
their dependencies as a [Docker
image](https://github.com/formbio/laava/pkgs/container/laava).

To get started, create a JSON file with your parameter values, similar to
[params-local-ss.json](https://raw.githubusercontent.com/formbio/laava/main/params-local-ss.json)
in this repo, and run it with:

```
nextflow run -profile local -params-file <your-params-file.json> main.nf
```


### Interactive Docker (`laava` container image)

For exploratory analysis or troubleshooting, you can also run the `laava` docker image
directly on the command line as an interactive container.

Assuming you have Docker installed, fetch the container image:

```
docker pull ghcr.io/formbio/laava:latest
```

Then run it interactively in your current directory:

```
docker run -v $(pwd):/data -w /data -it ghcr.io/formbio/laava:latest bash
```

### Python Runner

LAAVA now includes a unified Python-based runner that orchestrates the entire analysis pipeline:

```bash
python src/laava_run.py --input-bam <input.bam> --vector-fasta <vector.fasta> \
    --annotation-bed <annotation.bed> --itr-label ITR
```

This runner coordinates vector type analysis, read mapping, and alignment analysis in a single command. A minimal conda environment (laava_minimal.conda_env.yml) is also available for running the core analysis components.

Example using the test data:
```bash
# Using the included test samples
python src/laava_run.py \
    --input-bam test/samples/ss.subsample005.bam \
    --vector-fasta test/samples/ss.construct.fasta \
    --annotation-bed test/samples/ss.annotation.bed \
    --itr-label ITR \
    --output-dir test_runner_results/ss \
    --sample-id ss
```

## For developers

You can directly download or clone the repo to use the scripts directly.

```
$ git clone https://github.com/formbio/laava.git
```

There are several ways to satisfy the script dependencies locally.


### Option 1: Development docker image (`laava_dev` container image)

The `laava_dev.dockerfile` in this repo installs the scripts' dependencies, but not the
scripts themselves, into a Docker container image that you can then use to run the local
copies of the scripts. This allows you to edit the code in this repo in-place and run it
within the container environment without rebuilding the container.

To build the container image with the name `laava_dev` (you can use another name if you prefer):

```
docker build -t laava_dev:latest -f laava_dev.dockerfile .
```

To run the container in the current working directory:

```
docker run -v $(pwd):$(pwd) -w $(pwd) -it laava_dev:latest bash
```

This opens a Bash shell with the scripts in the PATH, and the original working directory mounted in place.


### Option 2: Conda


The conda (or mamba) channels and dependencies are in the configuration files
`laava.conda_env.yml`, `laava_dev.conda_env.yml`, and `laava_minimal.conda_env.yml`. 
The first two environments are similar and both will work for running LAAVA itself 
(with `_dev` including additional developer tools), while `_minimal` provides just 
the core components needed for the Python runner.

The conda channels and dependencies are in the configuration file `conda_env.yml`.
With this environment and a LaTeX installation (via e.g. apt), you'll have all the
dependencies you need to run LAAVA scripts directly on Linux, and *nearly* everything
you need on Mac.
The same environment is also used by the Docker container images internally.


First, install conda via [Miniconda or
Anaconda](https://www.anaconda.com/download/success).

Next, use the YAML configuration file to create a new conda environment and install its dependencies:

```
$ conda env create -f conda_env.yml
```

Finally, once installation completes, activate the new environment:

```
$ source activate laava
```

At this point the prompt should change to `(laava) $` and the executable scripts should
be available in your PATH.


## Testing

This repo includes small test files based on [PacBio's public AAV sequencing
examples](https://downloads.pacbcloud.com/public/dataset/AAV/).

### Automated local tests

The `test/` subdirectory in this repo contains small example input files and a Makefile
to run the scripts to reanalyze them and produce example HTML and PDF reports.

If you have Docker, Nextflow, and Make available, you can run a variety of tests from
the top directory of this repo.

* `make test` -- run both test samples using the Docker image directly, skipping Nextflow, and check the results quantitatively.
* `make sc` -- run the example self-complementary AAV (scAAV) sample with the Nextflow pipeline. This takes about 1-2 minutes.
* `make ss` -- run the example single-stranded AAV (ssAAV) sample. This takes about 2-3 minutes, including an additional flip/flop analysis step.
* `make all` -- run both example AAV samples using Nextflow.
* `make min` -- run the scAAV sample with the minimum number of required parameters, exercising the default behavior including guessing the construct vector type (sc/ss).
* `make folder` -- run both samples via folder input, exercising a batch processing mode.

Each of these commands will generate example HTML and PDF reports from the test datasets
included in the repo, which you can view locally.

