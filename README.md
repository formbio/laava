# LAAVA: Long-read AAV Analysis

A software package and executable bioinformatics workflow for the analysis of recombinant adeno-associated virus (rAAV) products by PacBio long-read sequencing.

For methods and interpretation, see: https://github.com/formbio/laava/wiki/Design-and-definitions

## Installation & Usage

Please read the [AAV tutorial](https://github.com/Magdoll/AAV/wiki/Tutorial:-Analyzing-AAV-Data)

### Standard Nextflow (recommended)

This code can be run as a standard [Nextflow](https://www.nextflow.io/) workflow.

When run this way, the workflow will automatically pull in the analysis scripts and
their dependencies as a [Docker
image](https://github.com/formbio/laava/pkgs/container/laava).

To get started, create a JSON file with your parameter values, similar to
[params-local-small.json](https://raw.githubusercontent.com/formbio/laava/main/params-local-small.json)
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
docker run -v $(pwd):$(pwd) -w $(pwd) -it ghcr.io/formbio/laava_dev:latest bash
```


## For developers

You can directly download or clone the repo to use the scripts directly.

```
$ git clone https://github.com/formbio/laava.git
```

There are several ways to satisfy the script dependencies locally.


### Option 1: Development docker image (`laava_dev` container image)

The `laava_dev.dockerfile` in this repo installs the scripts' dependencies, but not the
scripts themselves, into a Docker container image that you can then use to run the
local copies of the scripts.

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
`laava.conda_env.yml` and `laava_dev.conda_env.yml`. These two environments are similar,
and both will work for running LAAVA itself, but `_dev` includes some additional
developer tools.

First, install conda via Miniconda or Anaconda. Then, for example, suppose you have
[anaconda](https://docs.anaconda.com/anaconda/install/linux/) installed and the binary
is in `$HOME/anaCogentPy37/bin`. To make the installed scripts available in your
environment, you would add the binary to $PATH if it isn't there already:

```
$ export PATH=$HOME/anaCogentPy37/bin:$PATH
```

Next, use the YAML configuration file to create a new conda environment and install its dependencies:

```
$ conda env create -f laava.conda_env.yml
```

Finally, once installation completes, activate the new environment:

```
$ source activate laava
```

At this point the prompt should change to `(laava) $` and the executable scripts should be available in your PATH.


#### Option 3: Manual installation from scratch

Bypassing all the above, you can use other package managers to install the dependencies
individually.

The prerequisites to run these scripts include:

* Python 3.7 or later
* R 3.6 or later

Python packages:
* [biopython](https://anaconda.org/bioconda/biopython)
* [pysam](https://anaconda.org/bioconda/pysam)
* [parasail-python](https://anaconda.org/bioconda/parasail-python)

R packages:
* tidyverse
* flextable
* Rmarkdown


## Testing

The `test/` subdirectory in this repo contains a Makefile that can fetch example PacBio
datasets from a public server and run the scripts on them to reanalyze them and produce
example HTML and PDF reports.

Once you've completed installation (above), activate your conda environment or Docker container and change to the test directory:

```
cd test
```

To generate the HTML and PDF reports from the test dataset included in the repo (this takes about 1 minute):

```
make
```
