#### python modules ####

import yaml
import sys
import itertools


#### LOAD CONFIG FILE AND BARCODES YAML ####
# load config file specifiying samples and parameters
configfile: "config.yml"

# construct a dictionary of samples and load the yaml for each
sample_barcodes = {}
filenames = {}
for sample in config:
	#read barcode yaml
	with open(config[sample]["barcodes"], 'r') as handle:
		try:
			#load yaml
			sample_barcodes[sample] = yaml.safe_load(handle)
			
		# handle any yaml errors
		except yaml.YAMLError as exc:
			print(exc)
			
	#check that number of starts specified in config matches number of sets of barcodes
	#specified in barcode yaml
	if not(len(sample_barcodes[sample]) == len(config[sample]["barcodes_start"])):
		raise InputError(f"Number of sets of barcodes in yaml must match number of starts specified in config file for sample {sample}")
	
	#populate filenames dict with filenames for each barcode set
	filenames[sample] = {}
	filenames[sample]["counter"] = 0 #a counter for later on
	current_filenames = [""] #to store  filenames for the current set of barcodes
	
	#add to filenames for each set of barcodes
	for set_name in config[sample]["barcodes_start"].keys():
	
		#add each barcode to each filename
		current_filenames = [string + "_" for string in current_filenames]
		current_filenames = list(itertools.product(current_filenames, sample_barcodes[sample][set_name]))
		current_filenames = [ "".join(list) for list in current_filenames ]
				
		#add to filenames dict
		filenames[sample][set_name] = current_filenames 
		
		
#### target files ####
rule all:
	input: 
		expand("out/{sample}/counts.txt", zip, sample = config.keys())
		
#### merging ####

#merge overlapping reads
rule merge:
	input:
		r1 = lambda wildcards: config[wildcards.sample]["path"] + wildcards.sample + config[wildcards.sample]["R1_pattern"],
		r2 = lambda wildcards: config[wildcards.sample]["path"] + wildcards.sample + config[wildcards.sample]["R2_pattern"]
	output:
		merged = "out/{sample}/{sample}.merged.fastq.gz",
		proc_r1 = "out/{sample}/{sample}.unmerged_R1.fastq.gz",
		proc_r2 = "out/{sample}/{sample}.unmerged_R2.fastq.gz"
	params:
		A = lambda wildcards: config[wildcards.sample]["adapter1"],
		B = lambda wildcards: config[wildcards.sample]["adapter2"]
	shell:
		"""
		echo {params.A}
		echo {params.B}
		bbmerge.sh in1="{input.r1}" in2="{input.r2}" out="{output.merged}" outu1="{output.proc_r1}" outu2="{output.proc_r2}" adapter1={params.A} adapter2={params.B}
		"""

#retain only reads that are the correct length
rule filter:
	input:
		"out/{sample}/{sample}.merged.fastq.gz"
	output:	
		"out/{sample}/{sample}.merged.filtered.fastq.gz"
	params:
		len = lambda wildcards: config[wildcards.sample]["amplicon_length"]
	shell:
		"""
		bbduk.sh in={input} out={output} minlen={params.len} maxlen={params.len}
		"""

#run script to count barcodes
rule count:
	input:
		"out/{sample}/{sample}.merged.filtered.fastq.gz"
	output:
		"out/{sample}/{sample}.counts.txt"
	params:
		barcodes = lambda wildcards: sample_barcodes{wildcards.sample},
		starts = lambda wildcards: config[wildcards.sample]["capsid_start"],
		unzipped_reads = lambda wildcards, input: str(input)[:-3],
	shell:
		"""
		gunzip -f {input}
		perl src/barcodes.pl --barcodes {params.bc_file} --reads {params.unzipped_reads} --out {output} --start {params.start} --length {params.length} --fwdPrimer {params.prim}
		"""