# script to count barcodes in NGS data
# a barcode is a sequence within a read that we want to count instances of
# there may be more than one set barcodes per read
#
# barcodes can be 'variable' or 'constant'
# a 'constant' barcode is specified by a list of sequences (a 'set')
# matching can either be exact or allow mismatches
# a 'variable' barcode is specified by constant sequences either side of the barcode
# the set of barcodes are of variable length and sequence (eg insertion libraries)
#
# each set of barcodes is specified by a yaml file, and its start coordinate within the read must be given
# output file with counts for each barcode by set and combination of barcodes
# allow for mismatches and variability in barcode location

import re
import argparse
import sys
import yaml

# for reverse complmenting
tab = str.maketrans("ACTG", "TGAC")

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='Count barcodes in NGS reads')
	parser.add_argument('--fastq', '-f', help='Fastq file containing reads',  required=True)
	parser.add_argument('--barcodes', '-b', help='Barcodes yaml file specifying barcodes to find', required=True)
	parser.add_argument('--fPrimer', '-p', help='Constant region of forward primer used for PCR (must be common to all reads)', type=str, required=False)
	parser.add_argument('--out', '-o', help='Output file')
	parser.add_argument('--mismatches', '-m', help='Number of mismatches to allow')
	args = parser.parse_args()

	# parse barcodes yaml
	barcs = parse_barcs_yaml(args)
	print("barcodes yaml:")
	print(barcs)
	print("")
	
	# construct search strategy
	search = construct_search(barcs, args)
	print("search strategy:")
	print(search)
	print("")
	
	# count barcodes
	print("found barcodes:")
	counts = count_barcodes(args, search)

def count_barcodes(args, search):
	"""
	count barcodes that are specified barcs within reads in fastq file specified in args.fastq
	"""	
	check_count = 0
	#open fastq file and read every second line of four (lines with sequences)
	with open(args.fastq) as handle:
		for i, line in enumerate(handle):
			if i % 4 == 1:
				if args.fPrimer is not None:
				# if forward primer is specified, check for it in read
					dropped_count = 0
					if re.search(args.fPrimer, line):
						pass
					else:
					# if we can't find it, reverse-complement and check again
						reversed = reverse_complement(line)
						if re.search(args.fPrimer, reversed):
							print(f'reversed line {i}')
							line = reversed
						# if we still can't find it, skip this line
						else:
							dropped_count += 1
							continue
				# check for barcodes
				print('searching line:')
				found_barcs = find_barcodes_in_line(line, search)
				print(found_barcs)
				check_count += 1
			
				
			
			if check_count > 20:
				exit()
			
	if args.fPrimer is not None:
		print(f"dropped {dropped_count} read(s) because forward primer could not be identified in forward or reverse orientation")			
	
def construct_search(barcodes, args):
	"""
	Construct a list that specifies how to search for barcodes
	List contains one entry for each set of barcodes
	Each entry in the list is a dictionary, which has one key:value pair specifying the type of search, and a second containing information needed to make the search (a 'search dict')
	
	If the barcodes are 'variable', the 'type' is 'variable', and the search dict contains only one key:value pair:
		the key is 'regex' and the value is the regex with capture group to use
		
	If the barcodes are 'constant', and no mismatches are allowed, the 'type' is 'constant_exact', and the search dict
		contains the barcodes where the barcode names are the keys and the barcode sequences are the values
		
	If the barcodes are 'constant' and mistmatches are allowed, the 'type' is 'constant_regex' and the search dict
		contains the barcode names as keys and regexes used to search for the barcodes as values
	 
	"""
	search = []
	for i, set in enumerate(barcodes):
		name = list(barcodes[i].keys())[0]
		# if type is variable, construct regex to match 
		if barcodes[i][name]['type'] == 'variable':
			search_dict = {'type':'variable'}
			search_dict['forward'] = f"{barcodes[i][name]['before']}(.+){barcodes[i][name]['after']}"
			if args.fPrimer is None:
				search_dict['reverse'] = f"{reverse_complement(barcodes[i][name]['after'])}(.+){reverse_complement(barcodes[i][name]['before'])}"
			search.append(search_dict)
		# if type is constant, we need to check if we are allowing mismatches or not
		elif barcodes[i][name]['type'] == 'constant':
			# if number of mismatches is specified
			search.append(create_barcodes_search_dict(barcodes[i][name], args))
						
	return search

def create_barcodes_search_dict(barcodes_dict, args):
	"""
	create a dictionary of barcodes for constructing a search
	if mismatches == 0, this is just a dict where each key and value is one barcode
	if mismatches > 0, this is a dict where each key is a barcode and each value is a regex used to search for that barcode
	"""	
	search_dict = {}
	
	#get dictionary with names:seq for barcodes
	forward_barcs = barcodes_dict['barcodes']
	
	# get expected start and stop for barcodes
	search_dict['start'] = barcodes_dict['start']
	search_dict['stop'] = search_dict['start'] + len(list(forward_barcs.values())[0])
	
	# get number of mismatches
	if 'mismatches' in barcodes_dict:
		mismatches = barcodes_dict['mismatches']
	else:
		mismatches = 0
	
	# construct dictionary for regexes/exact matches
	if mismatches == 0:
		# if not allowing mismatches, just use sequences from yaml file
		search_dict['type'] = 'constant_exact'
		search_dict['forward_search'] = {key:value for key, value in forward_barcs.items()}
		if args.fPrimer is None:
			search_dict['reverse_search'] = {key:reverse_complement(value) for key, value in forward_barcs.items()}
	if mismatches > 0:
		# if allowing mismatches, use regexes
		search_dict['type'] = 'constant_regex'
		search_dict['forward_search'] = {key:create_mismatches_regex([value], mismatches) for key, value in forward_barcs.items()}
		if args.fPrimer is None:
			search_dict['reverse_search'] = {key:create_mismatches_regex([reverse_complement(value)], mismatches) for key, value in forward_barcs.items()}
		
	return search_dict
	
def create_mismatches_regex(sequence_list, mismatches):
	"""
	Create regex consisting that will match the 'sequence' with 'mismatches' number of mismatches
	""" 
	# if we've already done all the mismatches, return the list
	if mismatches == 0:
		return "|".join(sequence_list)
	# otherwise, take every element in the list, and add to the list new strings in which '.' replaces each letter
	else:
		new_sequence_list = []
		for seq in sequence_list:
			for i in range(len(seq)):
				new_sequence_list.append(seq[:i] + '.' + seq[i + 1:])
		# then call the function with one fewer mismatches
		return create_mismatches_regex(new_sequence_list, mismatches - 1)
		
def find_barcodes_in_line(line, search):
	"""
	look for the barcodes specified in 'search' in the sequence 'line'
	return a list of the names of the barcodes found
	"""
	# iterate over sets to search for
	found_barcodes = []
	for set in search:
		if set['type'] == 'constant_exact':
			matches = []
			# get part of read to check
			start = set['start']
			stop = set['stop']
			subread_forward = line[start:stop]
			
			# check if the subread matches any of the barcodes
			for name, barcode in set['forward_search'].items():
				if barcode == subread_forward:
					matches.append(name)
			
			# if we don't know what the forward primer is, we also need to check the reverse part of the read
			if 'reverse_search' in set:
				subread_reverse = reverse_complement(line[-start:-stop])
				for name, barcode in set['reverse_search'].items():
					if barcode == subread_reverse:
						matches.append(name)
			
			# check how many matches we found
			if len(matches) == 0:
				found_barcodes.append('none')
			elif len(matches) == 1:
				found_barcodes.append(matches[0])
			elif len(matches) > 1:
				found_barcodes.append('ambiguous')
			
		elif set['type'] == 'constant_regex':
			matches = []
			
			# get part of read to check
			start = set['start']
			stop = set['stop']
			subread_forward = line[start:stop]
			
			# check if subread matches any of the barcodes
			for name, barcode in set['forward_search'].items():
				if re.findall(barcode, subread_forward):
					matches.append(name)
			
			# if we don't know what the forward primer is, we also need to check the reverse part of the read
			if 'reverse_search' in set:
				subread_reverse = reverse_complement(line[-start:-stop])
				for name, barcode in set['reverse_search'].items():
					if re.findall(barcode, subread_reverse):
						matches.append(name)
						
			# check how many matches we found
			if len(matches) == 0:
				found_barcodes.append('none')
			elif len(matches) == 1:
				found_barcodes.append(matches[0])
			elif len(matches) > 1:
				found_barcodes.append('ambiguous')
			
		elif set['type'] == 'variable':
			matches = []
			fwd_regex = set['forward']
			
			# match regex
			matches = re.findall(fwd_regex, line)
			
			# match in reverse orientiation if necesary
			if 'reverse' in set:
				matches = matches + re.findall(set['reverse'], reverse_complement(line))
			
			# check how many matches we found
			if len(matches) == 0:
				found_barcodes.append('none')
			elif len(matches) == 1:
				found_barcodes.append(matches[0])
			elif len(matches) > 1:
				found_barcodes.append('ambiguous')	
		
	return found_barcodes
		
	
def parse_barcs_yaml(args):
	"""
	parse barcodes yaml file and check that it makes sense
	"""
	
	# read barcodes yaml
	with open(args.barcodes, 'r') as stream:
		try:
			barcodes = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)
			
	print(f"found {len(barcodes)} sets of barcodes")
	names = []
	
	# check that barcodes are specified correctly
	for i in range(len(barcodes)):
		try:
			name = list(barcodes[i].keys())[0]
			names.append(name)
			# if this set is constant, check it contains a start and some sequences
			if barcodes[i][name]['type'] == 'constant':
				num_barcs = len(barcodes[i][name]['barcodes'])
				start = barcodes[i][name]['start']
				
				# check that all the barcodes are the same length
				barcodes_list = list(barcodes[i][name]['barcodes'].values())
				if not(all([len(barc) == len(barcodes_list[0]) for barc in barcodes_list])):
					raise ValueError(f"all barcodes in constant set {name} must be the same length")
					
				# check that if mismatches are specified, the number is not greater than the length of the barcode
				if 'mismatches' in barcodes[i][name]:
					if barcodes[i][name]['mismatches'] > len(barcodes_list[0]):
						raise ValueError(f"number of mismatches must not be greater than the length of barcodes in set {name}")
				
				print(f"set {name} contains {num_barcs} barcodes, and starts at position {start} in read")
				
			# if this set is variable, check it contains a before and after sequence
			elif barcodes[i][name]['type'] == 'variable':
				before = barcodes[i][name]['before']
				after = barcodes[i][name]['after']
				print(f"set {name} will consist of all the sequences occuring between sequences {before} and {after} in the read")
		except KeyError:
			print("check barcodes yaml file is valid:")
			print("'constant' barcode sets must contain a 'start' and a list of barcodes")
			print("'variable' barcode sets must specify the sequence 'before' and 'after' the variable sequence")
			raise ValueError("please specify a valid barcodes yaml")
	
	# check that there are no duplicate names		
	if len(names) != len(set(names)):
		raise ValueError("Please specify unique names for each barcode set")
		
	return barcodes
	
def reverse_complement(seq):
    return seq.translate(tab)[::-1]
	
if __name__ == "__main__":
	main(sys.argv[1:])