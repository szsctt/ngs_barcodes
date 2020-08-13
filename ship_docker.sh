#!/bin/bash


rm -rf barcodes_docker
mkdir -p "barcodes_docker/config"

cp run.sh run_sim_analysis_docker.sh barcodes_docker
cp config/analysis_docker.yml config/barcodes_*.yml config/sim_*.yml barcodes_docker/config

now=$(date +"%Y_%m_%d")

tar -czvf ${now}_barcodes_docker.tar.gz barcodes_docker

rm -r barcodes_docker