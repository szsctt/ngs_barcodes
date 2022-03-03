#!/bin/bash

set -euo pipefail

cp ../setup/env.yml .
cp ../src/barcodes.py .
cp ../src/sim.py .
cp ../Snakefile .
cp ../run.sh .

#### singularity build ####
# data must be in data/ (in local directory) - automatically bind-mounted with singularity
# output is to out/ (in local directory) - automatically bind-mounted with singularity
# snakefile is at /usr/src/Snakefile (inside container)
# scripts are in /usr/src/ (inside container)

perl -pe "s/config\[wildcards.sample\]\[\'path\'\]/\'data\/\'/g" Snakefile > Snakefile.bak
perl -pe "s/src/\/usr\/src/g" Snakefile.bak > Snakefile

# micromaba updates a env called 'base'
perl -pe "s/barcodes2/base/g" env.yml > env.yml.bak
mv  env.yml.bak  env.yml


docker build . -t szsctt/barcodes:latest_singularity -t szsctt/barcodes:6_singularity -t szsctt/barcodes:latest -t szsctt/barcodes:6

#docker push szsctt/barcodes:6_singularity
#docker push szsctt/barcodes:latest_singularity
#docker push szsctt/barcodes:latest

#### docker build ####
# data mirrored to /usr/local/src/data/
# output mirrored to /usr/local/src/out/
# snakefile and scripts same as singularity build

perl -pe "s/data\//\/usr\/local\/src\/data\//g" Snakefile > Snakefile.bak
perl -pe "s/out\//\/usr\/local\/src\/out\//g" Snakefile.bak > Snakefile


docker build . -t szsctt/barcodes:latest_docker -t szsctt/barcodes:6_docker 

#docker push szsctt/barcodes:6_docker
#docker push szsctt/barcodes:latest_docker
