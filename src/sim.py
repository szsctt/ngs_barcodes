# script to simulate barcodes in NGS data
# a barcode is a sequence within a read that we want to count instances of
# there may be more than one set barcodes per read
#
# barcodes can be 'variable' or 'constant'
# a 'constant' barcode is specified by a list of sequences (a 'set')
# matching can either be exact or allow mismatches
# a 'variable' barcode is specified by constant sequences either side of the barcode
# the set of barcodes are of variable length and sequence (eg insertion libraries)
#
# each set of barcodes is specified by a yaml file
# output file with counts for each barcode by set and combination of barcodes
# allow for mismatches and variability in barcode location

from sys import argv
import argparse
import yaml
import numpy as np
import csv
import pdb

# for reverse translation
old = "acgtACGT"
new = "tgcaTGCA"
trans_table = str.maketrans(old,new)

def main(argv):
	parser = argparse.ArgumentParser(description='Simulate barcodes in NGS reads')
	parser.add_argument('--barcodes', '-b', help='Barcodes yaml file specifying reads and barcodes to use', required=True)
	parser.add_argument('--out-fastq_1', '-f1', help='Output R1 fastq file', default="sim.R1.fastq")
	parser.add_argument('--out-fastq_2', '-f2', help='Output R2 fastq file', default="sim.R2.fastq")
	parser.add_argument('--out-info', '-i', help='Output fastq file', default="sim_info.txt")	
	parser.add_argument('--read-len', '-l', help='length of simulated reads', default=150, type=int)	
	parser.add_argument('--n-sim', '-n', help='number of reads to simulate', default=100, type=int)
	parser.add_argument('--seed', '-s', help='seed for random number generator', default=12345, type=int)
	args = parser.parse_args()
	
	# parse and check input yaml file
	barcs = parse_yaml(args.barcodes, args.read_len)
	
	# open output files and simulate reads
	sim_reads(barcs, args)
	 		
	
	
	
def check_constant(entry):
	"""
	check required fields for set of constant barcodes
	"""
	assert 'barcodes' in entry
	
	assert len(entry['barcodes']) > 0
	
	# check each barcode has a sequence
	prob_sum = 0
	names = []
	lengths = []
	seqs = []
	for barc in entry['barcodes']:
		# make sure there's a sequence
		assert 'seq' in entry['barcodes'][barc]
		
		# make sure sequences are unique
		assert entry['barcodes'][barc]['seq'] not in seqs
		seqs.append(entry['barcodes'][barc]['seq'])
		
		# make sure lengths of the barcodes are all the same
		if len(lengths) == 0:
			lengths.append(len(entry['barcodes'][barc]['seq']))
		else:
			assert all([len(entry['barcodes'][barc]['seq']) == i for i in lengths])
		
		#make sure names are unique
		assert barc not in names
		names.append(barc) 
		
		# make sure there's a probability, and keep track of the sum of the probabilities
		if 'prob' not in entry['barcodes'][barc]:
			barc['prob'] = 1/len(entry['barcodes'])
		prob_sum += entry['barcodes'][barc]['prob']
		
	# check sum of probabilities is 1
	assert prob_sum == 1
	
	# return length this set will add to amplicon
	return lengths[0] + len(entry['after'])
	
def check_variable(entry):
	"""
	check required fields for a variable set of barcodes
	"""
	
	assert 'max_len' in entry
	assert isinstance(entry['max_len'], int)
	assert entry['max_len'] >= 0
	
	if 'min_len' in entry:
		assert isinstance(entry['min_len'], int)
		assert entry['max_len'] >= entry['min_len']
		assert entry['min_len'] >= 0
		
	else:
		entry['min_len'] = 0
		
	if 'len_multiple' in entry:
		assert isinstance(entry['len_multiple'], int)
		assert entry['len_multiple'] >= 0
		assert entry['max_len'] % entry['len_multiple'] == 0
		assert entry['min_len'] % entry['len_multiple'] == 0
	else:
		entry['len_multiple'] = 1
		
		
	# return length this set will add to amplicon
	min_len = entry['min_len'] + len(entry['after'])
	max_len = entry['max_len'] + len(entry['after'])
	
	return min_len, max_len
	
def parse_yaml(barcodes_path, read_len):
	"""
	check that yaml file at specified path is correct, and return it as a python object
	"""
	
	# read barcodes yaml
	with open(barcodes_path, 'r') as stream:
		try:
			barcodes = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)
			
	# check the yaml contains a list with at least one entry
	assert len(barcodes) > 0
	assert isinstance(barcodes, list)
	
	# check each entry in the list is valid
	set_names = []
	# only 'variable' sets can follow 'variable', sets, because 'constant' sets rely on position within the read
	seen_variable = False
	min_amp_length = 0
	max_amp_length = 0
	for i, entry in enumerate(barcodes):
	
		# check entry is a dict with one key/value pair
		# and that set names are unique
		assert len(entry) == 1
		name = list(entry.keys())[0]
		assert name not in set_names
		set_names.append(name)
		
		# required fields
		assert 'type' in entry[name]
		assert 'after' in entry[name]
		
		if i == 0:
			assert "before" in entry[name]
			min_amp_length += len(entry[name]['before'])
			max_amp_length += len(entry[name]['before'])
		
		# check specific fields
		if entry[name]['type'] == "constant":
			if seen_variable is True:
				raise ValueError("Sets of constant barcodes cannot follow sets of variable barcodes")
			barc_len = check_constant(entry[name])
			
			# check how much length this barcode will add to the amplicon
			max_amp_length +=  barc_len
			min_amp_length +=  barc_len
			
			
		elif entry[name]['type'] == "variable":
			min_len, max_len = check_variable(entry[name])
			max_amp_length +=  max_len
			min_amp_length +=  min_len
			
			
			seen_variable = True
		else:
			raise ValueError(f"type field of of {entry[0]} must be 'constant' or 'variable' (not {entry[0]['type']})")
			
	
	# minimimum length can't be less than the length of the reads
	# otherwise we don't have enough sequence to fill them
	assert min_amp_length >= read_len
	
	# if minimum length is more than 2x the read length, print a warning
	# because we won't be able to merge reads
	if min_amp_length >= 2 * read_len:
		print(f"WARNING: minimum fragment length is more than twice the read length.  Reads will not be able to be merged reliably.")
			
	
	return(barcodes)
	
def sim_reads(barcs, args):
	"""
	simulate reads and write output to specified files
	"""
	
	rng = np.random.default_rng(args.seed)
	
	# open output files
	with open(args.out_info, 'w', newline='') as info_file:
	
		# make csv writer
		infowriter = csv.writer(info_file, delimiter="\t")
		
		# write header to info file
		header = make_header(barcs)
		infowriter.writerow(header)
		
		with open(args.out_fastq_1, 'w') as fastq_1, open(args.out_fastq_2, 'w') as fastq_2:
		
			for id in range(args.n_sim):
				# generate fragment sequence
				frag, info = make_frag(barcs, rng)
				
				
				# 50% chance of reverse complementing read
				if rng.choice([True, False]):
					revcomp = True
					frag = reverse_complement(frag)
				else:
					revcomp = False
				
				# output information
				infowriter.writerow([id, revcomp] + info)
				
				# generate reads
				read_1 = frag[:args.read_len]
				read_2 = reverse_complement(frag[-args.read_len:])
				
				# write to fastq
				fastq_1.write(f"@read_{id}\n{read_1}\n+\n{''.join(['A']*len(read_1))}\n")
				fastq_2.write(f"@read_{id}\n{read_2}\n+\n{''.join(['A']*len(read_2))}\n")		
	
def make_header(barcs):
	"""
	construct header for information file, depending on the requested barcodes
	"""
	
	# information used for each type
	constant_info = ['set_name', 'start_position', 'type', 'name', 'seq']
	variable_info = ['set_name', 'start_position', 'type', 'length', 'seq']
	header = ['id', 'rev_comp']
	
	# iterate over sets of barcodes, and append information
	for barc_set in barcs:
		name = list(barc_set.keys())[0]
		if barc_set[name]['type'] == "constant":
			header = header + constant_info
		else:
			header = header + variable_info
		
	return header
	
def make_frag(barcs, rng):
	"""
	randomly generate a fragment and return the read and information about the barcodes in the fragment
	"""

	frag = [value['before'] for value in barcs[0].values()]
	info = []
	
	for barcs_set in barcs:
		# append name of set to info
		info = info + [name for name in barcs_set.keys()]
		name = info[-1]
		
		# add start_position to info
		info.append(sum([len(seq) for seq in frag]))
		
		# get a random constant barcode
		if barcs_set[name]['type'] == "constant":
			# append information
			info.append('constant')
		
			# get lists of names, sequences and probabilities
			names = [name for name in barcs_set[name]['barcodes'].keys()]
			seqs = [barcs_set[name]['barcodes'][i]['seq'] for i in names]
			probs = [barcs_set[name]['barcodes'][i]['prob'] for i in names]
			
			# get a random name
			barc_name = str(rng.choice(names, p = probs))
			info.append(barc_name)
			
			# append barcode sequence to frags
			frag.append(barcs_set[name]['barcodes'][barc_name]['seq'])
			info.append(barcs_set[name]['barcodes'][barc_name]['seq'])
			
			
		elif barcs_set[name]['type'] == "variable":
			# append information
			info.append('variable')
			
			# get available choices for number of bases
			lengths = range(barcs_set[name]['min_len'], barcs_set[name]['max_len'] + 1)
			lengths = [length for length in lengths if length % barcs_set[name]['len_multiple'] == 0]
			
			# get length
			length = rng.choice(lengths)
			info.append(length)
			
			# get random bases
			seq = "".join(np.random.choice(["A", "C", "G", "T"], length))
			info.append(seq)
			frag.append(seq)
		
		# append after sequence
		frag.append(barcs_set[name]['after'])
	
	
	return "".join(frag), info
	
def reverse_complement(seq):
	"""
	do reverse complement
	"""
	
	for base in seq:
		assert base in old
	
	return seq.translate(trans_table)[::-1]
	

if __name__ == "__main__":
	main(argv[1:])
	