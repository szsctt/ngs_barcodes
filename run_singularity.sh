#!/bin/bash
set -e

# script to run barcodes analysis pipeline
# reads must be stored in a folder called data/

# provide path to config file as command line argument
# this must be somewhere in the current directory, or a subdirectory of the current directory


# change the number below to set the max number of jobs run in parallel
JOBS=1

# usage: ./run.sh <config_file>

if [ "$#" -ne 1 ]; then
    echo "usage: ./run.sh <config_file>"
    exit 2
fi

# this is the location within the container that we will find the contents of the current directory
MOUNTDIR="/usr/local/src"
TMP="tmp"
mkdir -p $TMP
CONTAINER="szsctt/barcodes:v0.1.1"

# get name of container from string
TAG=$(echo $CONTAINER | cut -d':' -f2)
BASE=$(echo $CONTAINER | cut -d':' -f1)
CONTAINERNAME=$(echo $BASE | cut -d'/' -f2)
SIF="${CONTAINERNAME}_${TAG}.sif"

if [ ! -e "${SIF}" ]; then
    echo "pulling container..."
    singularity pull ${SIF} docker://${CONTAINER}
fi

# run container
cd $TMP && ln -sf ../$SIF ../config ../data .
SMK="snakemake --configfile ${MOUNTDIR}/${1} --jobs ${JOBS}"
singularity exec \
--bind "$(realpath ..):${MOUNTDIR}" \
--bind "$(pwd):/app/.snakemake" \
${SIF} \
bash -c "export PATH="/opt/conda/bin:$PATH" && cp -r /app/Snakefile /app/src . && ${SMK}"
cp -r out .. && cd .. && rm -r $TMP
