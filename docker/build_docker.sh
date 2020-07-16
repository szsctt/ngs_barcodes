#!/bin/bash

cp ../setup/env.yml .
cp ../Snakefile .
cp ../src/barcodes.py .

sudo docker build . -t szsctt/barcodes:latest -t szsctt/barcodes:1