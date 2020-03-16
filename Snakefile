
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
		temp("out/{sample}/{sample}.merged.filtered.fastq.gz")
	params:
		min_len = lambda wildcards: config[wildcards.sample]["min_length"],
		max_len = lambda wildcards: config[wildcards.sample]["max_length"]
	shell:
		"""
		bbduk.sh in={input} out={output} minlen={params.min_len} maxlen={params.max_len}
		"""
		
#### counting ####

#run script to count barcodes
rule count:
	input:
		"out/{sample}/{sample}.merged.filtered.fastq.gz"
	output:
		"out/{sample}_counts.txt"
	params:
		barcodes = lambda wildcards: config[wildcards.sample]["barcodes"],
		unzipped_reads = lambda wildcards, input: str(input)[:-3],
		prim = lambda wildcards: config[wildcards.sample]["fwdPrimer"] if "fwdPrimer" in config[wildcards.sample] else "",
		extend = lambda wildcards: config[wildcards.sample]["extend"] if "extend" in config[wildcards.sample] else 0,
		mismatches = lambda wildcards: config[wildcards.sample]["mismatches"] if "mismatches" in config[wildcards.sample] else 0
			
	shell:
		"""
		gunzip -f {input}
		perl src/barcodes.pl --barcodes {params.barcodes} --reads {params.unzipped_reads} --outpath out/{wildcards.sample}/ --mismatches {params.mismatches} --extend_search {params.extend} --prefix {wildcards.sample} --fwdPrimer {params.prim}
		(head -n 1 out/{wildcards.sample}/{wildcards.sample}_counts.txt && tail -n +2 out/{wildcards.sample}/{wildcards.sample}_counts.txt | sort) > out/{wildcards.sample}_counts.txt
		"""
		