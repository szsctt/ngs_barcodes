
#### load config file ####
# load config file specifiying samples and parameters
configfile: "config.yml"
				
#### target files ####
rule all:
	input: 
		expand("out/{sample}_counts.txt", sample = config.keys())
		
#### merging ####

#merge R1 and R2
rule merge:
	input:
		r1 = lambda wildcards: config[wildcards.sample]["path"] + wildcards.sample + config[wildcards.sample]["R1_pattern"],
		r2 = lambda wildcards: config[wildcards.sample]["path"] + wildcards.sample + config[wildcards.sample]["R2_pattern"]
	output:
		merged = temp("out/{sample}/{sample}.merged.fastq.gz"),
		proc_r1 = temp("out/{sample}/{sample}.unmerged_R1.fastq.gz"),
		proc_r2 = temp("out/{sample}/{sample}.unmerged_R2.fastq.gz")
	params:
		A = lambda wildcards: config[wildcards.sample]["adapter1"],
		B = lambda wildcards: config[wildcards.sample]["adapter2"]
	shell:
		"""
		bbmerge.sh in1="{input.r1}" in2="{input.r2}" out="{output.merged}" outu1="{output.proc_r1}" outu2="{output.proc_r2}" adapter1={params.A} adapter2={params.B}
		"""

#### filtering ####

#retain only reads that are the correct length
rule filter:
	input:
		"out/{sample}/{sample}.merged.fastq.gz"
	output:	
		temp("out/{sample}/{sample}.merged.filtered.fastq")
	params:
		min_len = lambda wildcards: config[wildcards.sample]["min_length"],
		max_len = lambda wildcards: config[wildcards.sample]["max_length"]
	shell:
		"""
		bbduk.sh in={input} out={output} minlen={params.min_len} maxlen={params.max_len}
		"""

#### counting ####

def translate_flag(wildcards):
	if "translate_insertion" in config[wildcards.sample]:
		if config[wildcards.sample]["translate_insertion"] is True:
			return "-t"
		else:
			return ""
	else:
		return ""

#run script to count barcodes
rule count:
	input:
		reads = "out/{sample}/{sample}.merged.filtered.fastq",
		barcodes = lambda wildcards: config[wildcards.sample]["barcodes"]
	output:
		"out/{sample}_counts.txt"
	params:
		barcodes = lambda wildcards: config[wildcards.sample]["barcodes"],
		prim = lambda wildcards: f"--fPrimer {config[wildcards.sample]['fwdPrimer']}" if "fwdPrimer" in config[wildcards.sample] else "",
		translate =  translate_flag
			
	shell:
		"""
		python3 src/barcodes.py --barcodes {params.barcodes} --fastq {input.reads} --out {output} {params.prim} {params.translate}
		"""
		