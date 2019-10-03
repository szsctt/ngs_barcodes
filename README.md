# Barcode counting

This is a pipeline for detection of barcodes in capsid sequences (Illumina).

It assumes data is paired-end, amplicon sequencing.  

## Installation

The pipeline is intended to be run on Mac OSX.  It may work on Linux, but has not been tested.  It may not work on Windows.

The pipeline requires `conda` to be installed.  Install either [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://docs.anaconda.com/anaconda/install/).

Once conda is installed, navigate to the installation directory and execute:

```
cd <install_directory>
setup/install.sh
```

This will check your conda installation and download all required tools.

## Inputs

Specify the inputs to the pipeline in a yaml file (an example `config.yml` is provided).
This config file should consist of nested dictionaries specifying where the data is and how it is to be processed.

### Sample

The config file must have one or more keys, each of which specify one set of reads.  The sample is the name of the file, without the suffixes (eg `_R1.fastq.gz` and `_R2.fastq.gz`).

### Reads

Provide a path to the folder containing the reads from this dataset (`path`).  Please do not use spaces or special characters in this (or any other) paths.

Also provide the suffix for each file (`R1_pattern`, `R2_pattern` eg `_R1.fastq.gz`, `_R2.fastq.gz`), the adapters used for sequencing (`adapter1` and `adapter2`), and the length of the amplicon (`amplicon_length`). If barcodes were added during amplification, note that the sequence of `fwdPrimer` must not contain those barcodes (in other words, it should be common to all the reads, regardless of which sample barcode was added).

### Barcodes

Sequences may contain any number of barcodes.  Provide barcodes for each sample in a seperate yaml (see example `example_barcodes.yml`).  This yaml is a nested dictionary with each outer key specifying one set of barcodes.  Each dictionary must contain a `start` key specifying the location in the read where that barcode is expected to start, and a number of key-value pairs with the name of each barcode (key) and its sequence (value).

### Optional parameters

There are a number of optional parameters which may be set for each sample.  These are `extend_search`, `mismatches` and `fwdPrimer`.  If a number of `mismatches` is specified, that number of mismatches will be allowed during the search. If the forward primer used for amplification (`fwdPrimer`) is specified (case insensitive), it will be used during barcode counting to check the orientation of the read.  A number is given for `extend_search`, the region of the read searched for the barcode will be extended that number of bases either side of the expected position of the barcode. 

## Running the pipeline

Once the config file and barcode yaml files has been correctly specified, run the pipeline from the installation directory:

```
cd <install_directory>
./barcodes
```

The results for each dataset are saved in sub-directories of the `out` directory.  For each sample, a merged and filtered fastq are saved, as well as counts for each barcode by set as well as all combinations of barcodes