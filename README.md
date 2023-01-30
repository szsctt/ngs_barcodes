# Barcode counting

This is a pipeline for detection of barcodes in capsid sequences (Illumina).

It assumes data is paired-end, amplicon sequencing.  

It can be run as a webapp, in a `docker` container, or using `snakemake`.

## Quickstart

To run an example using the docker container:

```
./run.sh config/
```

## Running the webapp

Navigate to the `ngs_barcodes` directory and run

```
docker compose up
```

Open a browser and go to `http://127.0.0.1:5000/`.


## Docker

Alternativley, run the pipeline in a docker container ([`szsctt/barcodes:latest`](https://hub.docker.com/r/szsctt/barcodes/tags)), using the `run.sh` script.  For this option, your data should be in a folder called `data`, and your config files should be a folder called `config`.  Your directory structure should look like this:

```
run.sh
data/
├─ my_reads_R1.fq.gz
├─ my_reads_R2.fq.gz
config/
├─ config.yml
├─ barcodes.yml
```

Run the pipeline using the `run.sh` script

```
./run.sh config/config.yml
```

## Using snakemake directly

The pipeline is intended to be run on Mac OSX.  It may work on Linux, but has not been tested.  It may not work on Windows.

The pipeline requires `conda` to be installed.  Install either [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://docs.anaconda.com/anaconda/install/).

Once conda is installed, navigate to the installation directory and execute:

```
cd <install_directory>
setup/install.sh
```

This will check your conda installation and download all required tools.


## Config file

Specify the inputs to the pipeline in a yaml file (an example `config.yml` is provided).
This config file should consist of nested dictionaries specifying where the data is and how it is to be processed. Data is assumed to be paired-end, amplicon Illumina sequencing.  There should be one 'block' in the yaml per sample

An example is:
```
sample1: # sample name
    path: "data/reads/" 
    R1_pattern: "_combined_R1.fastq"                 
    R2_pattern: "_combined_R2.fastq"                 
    adapter1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" 
    adapter2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 
    min_length: 152 
    max_length: 155     
    fwdPrimer: "ACCACCAGCACCAGAACCTGG"    
    translate_insertion: True     
    barcodes: "example_barcodes.yml"    
    
```


### Reads

Provide a path to the folder containing the reads from this dataset (`path`).  Please do not use spaces or special characters in this (or any other) paths.  The path can be the absolute path, or relative to the installation directory.

Also provide the suffix for each file (`R1_pattern`, `R2_pattern` eg `_R1.fastq.gz`, `_R2.fastq.gz`), the adapters used for sequencing (`adapter1` and `adapter2`), and the length of the amplicon (`amplicon_length`). If barcodes were added during amplification, note that the sequence of `fwdPrimer` must not contain those barcodes (in other words, it should be common to all the reads, regardless of which sample barcode was added).

### Barcodes

Sequences may contain any number of barcodes.  Provide barcodes for each sample in a seperate yaml (see example `example_barcodes.yml`).  This yaml is a nested dictionary with each outer key specifying one set of barcodes.  Each dictionary must contain a `start` key specifying the location in the read where that barcode is expected to start, and a number of key-value pairs with the name of each barcode (key) and its sequence (value).

### Forward primer

This is used to check each read is in the correct orientation.  It can be the common portion of the forward primer used for PCR prior to sequencing (i.e. remove any barcodes that are not common to all reads in the sample), or it can be any other sequence that is common to all the reads.

### Read length

In order to filter reads of only the expected lengths, specify `min_length` and `max_length`.

### translate\_insertion

If the sample contains variable barcodes (see below), specify if you would like this to be translated into amino acids in the output.  If the barcode length is not a multiple of three, it cannot be translated and will be output as a nucleotide sequence, enclosed in parentheses.

### Barcodes

Barcodes for each sample should be specified in a seperate yaml file. Each read may contain multiple sets of barcodes at different positions in the read.  Barcode sets can consist of either 'constant' or 'variable' barcodes

#### Constant barcodes

These are a part of the read that we would like to count instances of.  They are specified by sequence and location in the read.

An example yaml for one set of constant barcodes:
```
- set1: 	# name of this set of barcodes				
    type: constant
    start: 0 
    barcodes:
        barcode1: CGTAG
        barcode2: TCCTA
```

The barcode set name must begin the yaml block.  The (0-based) position where this barcode is expected within the read must also be specified (`start`).  The sequences of the barcodes (`barcodes`) should be specified, with each barcode given a name.  

One additional barcode name may be present in the output files: 'none', where none of the barcodes specified could be identified for a given read.

#### Variable barcodes

These are sequences that we want to count that are of variable sequence and length, but are flanked on either side by some constant sequence.  These are intended to address the problem of an insertion library in which many possible peptides, possibly of varying length, are inserted at some point in the read.

```
- set2:
    type: variable 
    before: GCCAATCA 
    after: GGAGCTTC 
```

Again, the set must have a name at the top, which is followed by the parameters of this set.  Instead of providing the sequence of each barcode, specify 'before' and 'after' sequences that immediately precede and follow the insertion site.  Any read in which the before and after sequences cannot be identified will be counted as 'none'. Any read in which the after sequence follows the before sequence (with nothing in between) will be counted as 'no\_insertion'.   

## Running the pipeline

Once the config file and barcode yaml files has been correctly specified, run the pipeline from the installation directory.  The number of cores to use must be specified.

```
cd <install_directory>
./barcodes --cores <cores>
```
Each job only uses one core, so using more tha one core will only speed the pipeline up if more than one sample is being run at once.

Snakemake requires that a number of cores be specified - increasing this will allow jobs to be run in parallel, but this will only speed up execution if more than one sample is being processed simultaneously.

Other snakemake options can also be passed in.  For example:
```
./barcodes --cores 1 --rerun-incomplete
```
will re-run any incomplete jobs from previous runs.

See the [snakemake docs](https://snakemake.readthedocs.io/en/stable/executable.html) for more information.

## Output

The results for each dataset are saved in sub-directories of the `out` directory.  For each sample, a merged and filtered fastq are saved, as well as the barcode counts (`out/{sample}_counts.txt`).

## Running the flask app

To run the flask app:

```
# start redis
redis-server

# start a redis worker
rq worker

# run flask app
flask run
```
