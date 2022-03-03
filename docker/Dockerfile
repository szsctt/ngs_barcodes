FROM mambaorg/micromamba:0.21.2

#  $ docker build . -t szsctt/barcodes:latest -t szsctt/barcodes:4
#  $ docker run --rm -it szsctt/barcodes:latest /bin/bash
#  $ docker push szsctt/barcodes:latest
#  $ docker push szsctt/barcodes:4


COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba install -y -f /tmp/env.yml &&\
	micromamba clean --all -y
 
COPY --chown=$MAMBA_USER:$MAMBA_USER Snakefile /usr/src
COPY --chown=$MAMBA_USER:$MAMBA_USER barcodes.py /usr/src/barcodes.py
COPY --chown=$MAMBA_USER:$MAMBA_USER sim.py /usr/src/sim.py
WORKDIR /usr/src
