from os import path

#### load config file ####
# load config file specifiying samples and parameters
# configfile: "config.yml" # supply on command line
	
	
#### wildcard constraints
wildcard_constraints:
	sample = "|".join(list(config.keys()))			
	
#### target files ####
rule all:
	input: 
		expand("out/{sample}_counts.txt", sample = config.keys()),
		"out/qc_input/multiqc_report.html",
		"out/qc_filt/multiqc_report.html"
		
rule fastqc_input:
	input:
		r1 = lambda wildcards: f"{path.normpath(config[wildcards.sample]['path'])}/{wildcards.sample + config[wildcards.sample]['R1_pattern']}",
		r2 = lambda wildcards: f"{path.normpath(config[wildcards.sample]['path'])}/{wildcards.sample + config[wildcards.sample]['R2_pattern']}"
	output:
		r1_zip =  temp("out/qc_input/{sample}_1_fastqc.zip"),
		r2_zip =  temp("out/qc_input/{sample}_2_fastqc.zip")
	params:
		# use temporary directories becuase running multiple samples in the same directory can set up race condition
		zip_path1=lambda wildcards, output: f"{path.dirname(output.r1_zip)}/{wildcards.sample}1/{wildcards.sample}{config[wildcards.sample]['R1_pattern'].split('.')[0]}_fastqc.zip",
		tempdir1 = lambda wildcards, output: f"{path.dirname(output.r1_zip)}/{wildcards.sample}1",
		zip_path2=lambda wildcards, output: f"{path.dirname(output.r2_zip)}/{wildcards.sample}2/{wildcards.sample}{config[wildcards.sample]['R2_pattern'].split('.')[0]}_fastqc.zip",
		tempdir2 = lambda wildcards, output: f"{path.dirname(output.r2_zip)}/{wildcards.sample}2",

	shell: 
		"""
		mkdir -p {params.tempdir1}
		mkdir -p {params.tempdir2}
		fastqc -o {params.tempdir1} {input.r1} 
		fastqc -o {params.tempdir2} {input.r2}
		mv {params.zip_path1} {output.r1_zip}
		mv {params.zip_path2} {output.r2_zip}
		rm -r {params.tempdir1} {params.tempdir2}
		"""
		
def multiqc_input(wildcards):
	files  = [f"out/qc_input/{samp}_1_fastqc.zip" for samp in config.keys()]
	files2 = [f"out/qc_input/{samp}_2_fastqc.zip" for samp in config.keys()]
	return files + files2
		
rule multiqc:
	input:
		multiqc_input
	output:
		"out/qc_input/multiqc_report.html"
	params:
		outdir = lambda wildcards, output: path.dirname(output[0])

	shell:
		"multiqc {input} -f --outdir {params.outdir}"
		
#### merging ####

#merge R1 and R2
rule merge:
	input:
		r1 = lambda wildcards: f"{path.normpath(config[wildcards.sample]['path'])}/{wildcards.sample + config[wildcards.sample]['R1_pattern']}",
		r2 = lambda wildcards: f"{path.normpath(config[wildcards.sample]['path'])}/{wildcards.sample + config[wildcards.sample]['R2_pattern']}"
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

rule fastqc_filtered:
	input:
		reads = "out/{sample}/{sample}.merged.filtered.fastq"
	output:
		reads_zip =  temp("out/qc_filt/{sample}.merged.filtered_fastqc.zip")
	params:
		# use temporary directories becuase running multiple samples in the same directory can set up race condition
		zip_path1=lambda wildcards, output: f"{path.dirname(output.reads_zip)}/{wildcards.sample}/{wildcards.sample}.merged.filtered_fastqc.zip",
		tempdir1 = lambda wildcards, output: f"{path.dirname(output.reads_zip)}/{wildcards.sample}",

	shell: 
		"""
		mkdir -p {params.tempdir1}
		fastqc -o {params.tempdir1} {input.reads} 
		mv {params.zip_path1} {output.reads_zip}
		rm -r {params.tempdir1}
		"""
		
rule multiqc_filtered:
	input:
			[f"out/qc_filt/{samp}.merged.filtered_fastqc.zip" for samp in config.keys()]
	output:
		"out/qc_filt/multiqc_report.html"
	params:
		outdir = lambda wildcards, output: path.dirname(output[0])

	shell:
		"multiqc {input} -f --outdir {params.outdir}"

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
		
