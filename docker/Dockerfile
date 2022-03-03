FROM mambaorg/micromamba:0.19.1

#  $ docker build . -t szsctt/barcodes:latest -t szsctt/barcodes:4
#  $ docker run --rm -it szsctt/barcodes:latest /bin/bash
#  $ docker push szsctt/barcodes:latest
#  $ docker push szsctt/barcodes:4


USER root

ADD env.yml /tmp/env.yml

RUN echo "test!"

RUN micromamba env update -n base -f /tmp/env.yml &&\
	micromamba clean --all -y
 
COPY --chown=$MAMBA_USER:$MAMBA_USER Snakefile /usr/src
COPY --chown=$MAMBA_USER:$MAMBA_USER barcodes.py /usr/src/barcodes.py
COPY --chown=$MAMBA_USER:$MAMBA_USER sim.py /usr/src/sim.py
WORKDIR /usr/src
