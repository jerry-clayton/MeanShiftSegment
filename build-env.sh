#!/usr/bin/env -S bash --login
##!/bin/bash

set -euo pipefail

basedir=$( cd "$(dirname "$0")" ; pwd -P)

echo "Creating conda environment from ${basedir}/ms-env.yml"
conda env create --name "ms-env" -f ${basedir}/ms-env.yml

source $(conda info --base)/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate ms-env

echo "Current Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

echo "Running  R package install script"
# Run the R script
Rscript --verbose ${basedir}/run/create-ms-env.R