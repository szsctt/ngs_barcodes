#  $ docker build . -t szsctt/barcodes:latest -t szsctt/barcodes:6
#  $ docker run --rm -it szsctt/barcodes:latest /bin/bash
#  $ docker push szsctt/barcodes:latest
#  $ docker push szsctt/barcodes:6

FROM mambaorg/micromamba:0.22.0
USER root

COPY setup/env.yml /tmp/env.yml

RUN micromamba update -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
 
# snakemake stuff
COPY Snakefile /app/Snakefile
COPY src /app/src

# flask app stuff
COPY setup.py /app/setup.py
COPY helpers.py /app/helpers.py
COPY app.py /app/app.py
COPY static /app/static
COPY templates /app/templates


WORKDIR /app

CMD ["flask", "run", "--host=0.0.0.0"]
