#!/bin/bash
set -e

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

# simulate some data
echo "simulating dataset 1, replicate 1..."
docker run --rm -v "$(pwd):${MOUNTDIR}" szsctt/barcodes:5_docker  \
python3 /usr/src/sim.py \
--barcodes "${MOUNTDIR}/config/sim_0.yml" \
--out-fastq_1 "${MOUNTDIR}/${DATADIR}/sim_0.R1.fq" \
--out-fastq_2 "${MOUNTDIR}/${DATADIR}/sim_0.R2.fq" \
--out-info "${MOUNTDIR}/${DATADIR}/sim_0.info.txt" \
--read-len 150 \
--n-sim 10000 \
--seed 12345

echo "simulating dataset 1, replicate 2..."
docker run --rm -v "$(pwd):${MOUNTDIR}" szsctt/barcodes:5_docker  \
python3 /usr/src/sim.py \
--barcodes "${MOUNTDIR}/config/sim_0.yml" \
--out-fastq_1 "${MOUNTDIR}/${DATADIR}/sim_1.R1.fq" \
--out-fastq_2 "${MOUNTDIR}/${DATADIR}/sim_1.R2.fq" \
--out-info "${MOUNTDIR}/${DATADIR}/sim_1.info.txt" \
--read-len 150 \
--n-sim 10000 \
--seed 23456

echo "simulating dataset 2, replicate 1..."
docker run --rm -v "$(pwd):${MOUNTDIR}" szsctt/barcodes:5_docker  \
python3 /usr/src/sim.py \
--barcodes "${MOUNTDIR}/config/sim_1.yml" \
--out-fastq_1 "${MOUNTDIR}/${DATADIR}/sim_2.R1.fq" \
--out-fastq_2 "${MOUNTDIR}/${DATADIR}/sim_2.R2.fq" \
--out-info "${MOUNTDIR}/${DATADIR}/sim_2.info.txt" \
--read-len 150 \
--n-sim 10000 \
--seed 34567

echo "simulating dataset 2, replicate 2..."
docker run --rm -v "$(pwd):${MOUNTDIR}" szsctt/barcodes:5_docker  \
python3 /usr/src/sim.py \
--barcodes "${MOUNTDIR}/config/sim_1.yml" \
--out-fastq_1 "${MOUNTDIR}/${DATADIR}/sim_3.R1.fq" \
--out-fastq_2 "${MOUNTDIR}/${DATADIR}/sim_3.R2.fq" \
--out-info "${MOUNTDIR}/${DATADIR}/sim_3.info.txt" \
--read-len 150 \
--n-sim 10000 \
--seed 45678


# analyse simulated data
docker run --rm -it -v "$(pwd):${MOUNTDIR}" szsctt/barcodes:5_docker \
snakemake \
--snakefile /usr/src/Snakefile \
--configfile "${MOUNTDIR}/config/analysis_docker.yml" \
--cores 2

