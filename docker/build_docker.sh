#!/bin/bash

cp ../setup/env.yml .
cp ../src/barcodes.py .
cp ../src/sim.py .


cp ../Snakefile .

#### singularity build ####
# data must be in data/ (in local directory) - automatically bind-mounted with singularity
# output is to out/ (in local directory) - automatically bind-mounted with singularity
# snakefile is at /usr/src/Snakefile (inside container)
# scripts are in /usr/src/ (inside container)

perl -pe "s/config\[wildcards.sample\]\[\'path\'\]/\'data\/\'/g" Snakefile > Snakefile.bak
perl -pe "s/src/\/usr\/src/g" Snakefile.bak > Snakefile


docker build . -t szsctt/barcodes:latest_singularity -t szsctt/barcodes:5_singularity -t szsctt/barcodes:latest -t szsctt/barcodes:5

docker push szsctt/barcodes:5_singularity
docker push szsctt/barcodes:latest_singularity
docker push szsctt/barcodes:latest

#### docker build ####
# data mirrored to /usr/local/data/
# output mirrored to /usr/local/out/
# snakefile and scripts same as singularity build

perl -pe "s/data\//\/usr\/local\/data\//g" Snakefile > Snakefile.bak
perl -pe "s/out\//\/usr\/local\/out\//g" Snakefile.bak > Snakefile


docker build . -t szsctt/barcodes:latest_docker -t szsctt/barcodes:5_docker 

docker push szsctt/barcodes:5_docker
docker push szsctt/barcodes:latest_docker
