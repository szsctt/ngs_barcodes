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


# pull docker image
echo "pulling docker image..."
docker pull szsctt/barcodes:5_docker

# data must always be stored in data/ (in folder called 'data' within current directory)
# for this example, we simulate some data and save it in /data
DATADIR="data/"
rm -rf ${DATADIR}
mkdir -p ${DATADIR}

# this is the location within the container that we will find the contents of the current directory
MOUNTDIR="/usr/local"

# analyse simulated data
echo 
echo "analysing data..."
docker run --rm -it -v "$(pwd):${MOUNTDIR}" szsctt/barcodes:5_docker \
snakemake \
--snakefile /usr/src/Snakefile \
--configfile "${MOUNTDIR}/${1}" \
--cores 2
