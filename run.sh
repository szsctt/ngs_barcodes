#!/bin/bash
set -e

# script to run barcodes analysis pipeline
# reads must be stored in a folder called data/

# provide path to config file as command line argument
# this must be somewhere in the current directory, or a subdirectory of the current directory

# usage: ./run.sh <config_file>

if [ "$#" -ne 1 ]; then
    echo "usage: ./run.sh <config_file>"
    exit 2
fi

# this is the location within the container that we will find the contents of the current directory
MOUNTDIR="/usr/local"

# analyse simulated data
echo 
echo "analysing data..."
docker run --rm -it -v "$(pwd):${MOUNTDIR}" szsctt/barcodes:5_docker \
snakemake \
--snakefile /usr/src/Snakefile \
--configfile "${MOUNTDIR}/${1}" \
--jobs 1

