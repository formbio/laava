#!/bin/bash
# Usage: run-workflow-local.sh
set -o pipefail
set -o errexit

SCRIPT_DIR=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# Run nextflow locally
nextflow run -profile local main.nf -params-file params-local-small-sc-no-ff.json
