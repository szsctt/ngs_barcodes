# Barcode counting

This is a pipeline for detection of barcodes in capsid sequences (Illumina).

It assumes data is paired-end, amplicon sequencing.  

## Installation

The pipeline is intended to be run on Mac OSX.  It may work on Linux, but has not been tested.  It is unlikely to work on Windows.

The pipeline requires `conda` to be installed.  Install either [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://docs.anaconda.com/anaconda/install/).

Once conda is installed, navigate to the installation directory and execute:
```
cd <install_directory>
setup/install.sh
```

This will check your conda installation and download all required tools.

## Inputs

Specify the inputs to the pipeline in a yaml file (an example `config.yml` is provided).
Each dataset may consist of paired-end reads from more than one sample.

### Reads

Provide a path to the folder containing the reads from this dataset (`path`).  Please do not use spaces or special characters in this (or any other) paths.

Also provide the suffix for each file (`R1_pattern`, `R2_pattern` eg `_R1.fastq.gz`, `_R2.fastq.gz`), the forward primer used for amplification (`fwdPrimer`), the adapters used for sequencing (`adapter1` and `adapter2`), and the length of the amplicon (`amplicon_length`).   If sample barcodes were added during amplification, note that the sequence of `fwdPrimer` must not contain those barcodes (in other words, it should be common to all the reads, regardless of which sample barcode was added).

### Barcodes

Sequences may contain two types of barcodes. Sample barcodes specify which sample a read originated from, and were added during amplification.  Capsid barcodes are unique to particular capsids. 

#### Sample barcodes
Optionally, if multiple samples were combined in one sequencing run, barcodes for each sample can be specified.  

In this case, specify the path to the file containing the barcodes and their names (`sample_barcodes`).  This file must have three columns: the first is a barcode ID (a number), the second is the barcode sequence (eg GATTAC) and the third is a barcode name.  Each column is separated by a space, and each line contains one barcode. 

#### Capsid barcodes
Provide a path to a file containing the barcodes and their names (`capsid_barcodes`).  This file must have three columns: the first is a barcode ID (a number), the second is the barcode sequence (eg GATTAC) and the third is a barcode name.  Each column is separated by a space, and each line contains one barcode.

Also provide the position within the amplicon where the barcode is located, by specifying the `capsid_start` and `capsid_length`.  Note that capsid\_start should specify the first base of the barcode after the removal of any sample barcodes.

## Running the pipeline

Once the config file has been correctly specified, run the pipeline from the installation directory:

```
cd <install_directory>
./barcodes
```

The results for each dataset are saved in sub-directories of the `out` directory.