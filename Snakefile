#### python modules ####

import glob
import os
import pathlib


#### CONFIG FILE ####
#which datasets are to be aligned to which hosts and viruses are specified in dsets.yaml
configfile: "config.yml"

#format of config file:
#dataset1:
#    path: "~/OneDrive - CSIRO/Projects/shuf/barcodes/test_data/"
#    capsid_barcodes: "~/OneDrive - CSIRO/Projects/shuf/barcodes/test_data/BC270619.txt"
#    sample_barcodes: None
#    R1_pattern: _R1.fastq.gz
#    R2_pattern: _R2.fastq.gz
#    adapter1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#    adapter2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

DATASETS = [] #to store datasets
SAMPLES = [] #to store samples

for dataset in config:
	#get files in directory and strip off specified suffix (eg "_1.fastq.gz")
	suffix = config[dataset]["R1_pattern"]
	p=pathlib.Path(config[dataset]["path"])
	print(f"searching for samples in {str(p)}/*{suffix}")
	samples = [os.path.basename(f)[:len(suffix)*-1] for f in p.glob(f"*{suffix}")]
	if not samples:
		print(f"can't find any samples for dataset {dataset}")
	for sample in samples:
		DATASETS.append(dataset)
		SAMPLES.append(sample)
		
		
#### target files ####
rule all:
	input: 
		expand("out/{dataset}/{sample}.counts.txt", zip, dataset = DATASETS, sample = SAMPLES)
		
#### merging ####

#merge overlapping reads
rule merge:
	input:
		r1 = lambda wildcards: config[wildcards.dataset]["path"] + wildcards.sample + config[wildcards.dataset]["R1_pattern"],
		r2 = lambda wildcards: config[wildcards.dataset]["path"] + wildcards.sample + config[wildcards.dataset]["R2_pattern"]
	output:
		merged = "out/{dataset}/{sample}.merged.fastq.gz",
		proc_r1 = "out/{dataset}/{sample}.unmerged_R1.fastq.gz",
		proc_r2 = "out/{dataset}/{sample}.unmerged_R2.fastq.gz"
	params:
		A = lambda wildcards: config[wildcards.dataset]["adapter1"],
		B = lambda wildcards: config[wildcards.dataset]["adapter2"]
	shell:
		"""
		echo {params.A}
		echo {params.B}
		bbmerge.sh in1="{input.r1}" in2="{input.r2}" out="{output.merged}" outu1="{output.proc_r1}" outu2="{output.proc_r2}" adapter1={params.A} adapter2={params.B}
		"""

#retain only reads 150 bp in length
rule filter:
	input:
		"out/{dataset}/{sample}.merged.fastq.gz"
	output:	
		"out/{dataset}/{sample}.merged.filtered.fastq.gz"
	shell:
		"""
		bbduk.sh in={input} out={output} minlen=150 maxlen=150
		"""
		
#run script to count barcodes
rule count:
	input:
		"out/{dataset}/{sample}.merged.filtered.fastq.gz"
	output:
		"out/{dataset}/{sample}.counts.txt"
	params:
		bc_file = lambda wildcards: config[wildcards.dataset]["capsid_barcodes"],
		start = lambda wildcards: config[wildcards.dataset]["capsid_start"],
		length = lambda wildcards: config[wildcards.dataset]["capsid_length"],
		prim = lambda wildcards: config[wildcards.dataset]["fwdPrimer"],
		unzipped_reads = lambda wildcards, input: str(input)[:-3],
	shell:
		"""
		gunzip -f {input}
		perl src/barcodes.pl --barcodes {params.bc_file} --reads {params.unzipped_reads} --out {output} --start {params.start} --length {params.length} --fwdPrimer {params.prim}
		"""