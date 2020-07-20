
#!/bin/bash

cp ../setup/env.yml .
cp ../src/barcodes.py .
cp ../src/sim.py .


cp ../Snakefile .
# be a bit more restrictive about where data and results are, to make sure
# that they're always in the directory in which the singularity sif file is
perl -pe "s/config\[wildcards.sample\]\[\'path\'\]/\'data\/\'/g" Snakefile > Snakefile.bak
perl -pe "s/src/\/src/g" Snakefile.bak > Snakefile



docker build . -t szsctt/barcodes:latest -t szsctt/barcodes:4
