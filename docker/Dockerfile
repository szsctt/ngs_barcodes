#  $ docker build . -t szsctt/barcodes:latest -t szsctt/barcodes:6
#  $ docker run --rm -it szsctt/barcodes:latest /bin/bash
#  $ docker push szsctt/barcodes:latest
#  $ docker push szsctt/barcodes:6

FROM mambaorg/micromamba:0.22.0
USER root

COPY env.yml /tmp/env.yml

RUN micromamba install -y -f /tmp/env.yml && \
    micromamba clean --all --yes
 
COPY Snakefile /usr/src/Snakefile
COPY barcodes.py /usr/src/barcodes.py
COPY sim.py /usr/src/sim.py
COPY run.sh /usr/src/run.sh

WORKDIR /usr/src

