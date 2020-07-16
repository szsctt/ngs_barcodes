FROM continuumio/miniconda:4.7.12

#  $ docker build . -t szsctt/barcodes:latest -t szsctt/barcodes:3
#  $ docker run --rm -it szsctt/barcodes:latest /bin/bash
#  $ docker push szsctt/barcodes:latest
#  $ docker push szsctt/barcodes:3

ADD env.yml /tmp/env.yml
RUN /opt/conda/bin/conda env update -n base -f /tmp/env.yml &&\
	/opt/conda/bin/conda clean --all -y 
	
ADD Snakefile .
ADD barcodes.py src/barcodes.py
