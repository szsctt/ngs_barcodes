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


import gzip
from mimetypes import guess_type
from functools import partial
import re
import argparse
import sys
from os import path
import yaml
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
from mimetypes import guess_type
from functools import partial
import csv
import pdb

# for reverse complmenting
old_chars = "ACTGactg"
new_chars = "TGACtgac"
tab = str.maketrans("ACTGactg", "TGACtgac")

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='Count barcodes in NGS reads')
	parser.add_argument('--fastq', '-f', help='Fastq file containing reads',  required=True)
	parser.add_argument('--barcodes', '-b', help='Barcodes yaml file specifying barcodes to find', required=True)
	parser.add_argument('--fPrimer', '-p', help='Constant region of forward primer used for PCR (must be common to all reads)', type=str, required=True)
	parser.add_argument('--out', '-o', help='Output file', default="counts.txt")
	parser.add_argument('--debug', help='Produce extra output useful for debugging', action='store_true')
	parser.add_argument('--debug-output', help="output directory for debugging output", default = "out/")
	args = parser.parse_args()

	# check arguments
	if args.debug is True:
		path.isdir(args.debug_output)
		# make csv

	# parse barcodes yaml
	barcs = parse_barcs_yaml(args)
	
	# construct search strategy
	search = construct_search(barcs, args)
	
	# count barcodes and write output
	if args.debug is False:
		counts = count_barcodes(args, search, False)	
	else:
		counts = count_barcodes(args, search, True, args.debug_output)
		
		
		# write info
	write_counts(args.out, counts, search)
	
	print(f"saved counts in file {args.out}")

def count_barcodes(args, search, debug=False, debug_read_folder = ""):
	"""
	count barcodes that are specified barcs within reads in fastq file specified in args.fastq
	"""	
	
	
	checked_count = 0
	dropped_count = 0
	rev_count = 0
	ambiguous_fPrimer = 0
	counts = {} # to store counts of combinations of barcodes
	
	if debug is True:
		info = {}
		debug_info = path.realpath(debug_read_folder + "/debug_info.tsv")
		debug_info_handle = open(debug_info, "w", newline = "")
		writer = csv.DictWriter(debug_info_handle, 
								fieldnames = ['read_name', 'dropped', 'reversed', 'barcodes'],
								delimiter = '\t')
		writer.writeheader()
		
		
	
	# open fastq file and read every second line of four (lines with sequences)
	# handle gzipped files as well as non-gzipped
	# https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz
	encoding = guess_type(args.fastq)[1]  # uses file extension
	_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

	with _open(args.fastq) as handle:
		for record in SeqIO.parse(handle, "fastq"):
			if debug is True:
				info['read_name'] = record.id
		
			# check for forward primer in read:
			reversed = False
			
			# check for forward primer in read - if not found, take reverse complement
			n_matches = record.seq.lower().count(args.fPrimer.lower())
			if n_matches == 0:
				# check reverse complement
				record.seq = record.seq.reverse_complement()
				n_matches = record.seq.lower().count(args.fPrimer.lower())
				reversed = True
				# if we still can't find it, drop this read
				if n_matches == 0:
					dropped_count += 1
					if debug is True:
						info['dropped'] = True
						info['reversed'] = 'NA'
						info['barcodes'] = 'NA'
					continue
					
			# check for multiple matches
			if n_matches > 1:
					ambiguous_fPrimer += 1
					dropped_count += 1
					continue
					
			if reversed is True:
				rev_count += 1
				if debug is True:
					info['reversed'] = True
					info['dropped'] = False
			else:
				if debug is True:
					info['reversed'] = False
					info['dropped'] = False
					
			# check for barcodes
			found_barcs = find_barcodes_in_line(str(record.seq), search)

			checked_count += 1
			if debug is True:
				info['barcodes'] = "__".join(found_barcs)
			
			# increment count for this combination
			if debug is True:
				writer.writerow(info)
				
			counts = increment_counter(counts, found_barcs)
	
	print(f"checked {checked_count} reads in total; {rev_count} of these were correctly reversed")
	print(f"dropped {dropped_count} read(s) because forward primer could not be identified in forward or reverse orientation")	
	print(f"forward primer appeared more than once in {ambiguous_fPrimer} reads: if this number is high, consider re-running with a longer forward primer sequence")	
	
	if debug is True:
		debug_info_handle.close()
		
	return counts
	
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
			
			# add type, name, and bool specifiying if we want to translate
			search_dict = {'type':'variable'}
			search_dict['name'] = name
			search_dict['trans'] =  barcodes[i][name]['translate']
			
			# if we allow mismatches
			if 'mismatches' in barcodes[i][name]:
				mismatches = barcodes[i][name]['mismatches']
			else:
				mismatches = 0
			
			# construct regex for searching
			search_dict['forward'] = construct_variable_regex(barcodes[i][name]['before'], barcodes[i][name]['after'], mismatches)
			
			#search_dict['forward'] = f"{barcodes[i][name]['before']}(.+){barcodes[i][name]['after']}"
			search.append(search_dict)
			
		# if type is constant, we need to check if we are allowing mismatches or not
		elif barcodes[i][name]['type'] == 'constant':
			# if number of mismatches is specified
			search_dict = create_barcodes_search_dict(barcodes[i][name], args)
			search_dict['name'] = name
			search.append(search_dict)
						
	return search
	
def construct_counts_access(barcs, dict_name):
	access = dict_name
	for barc in barcs:
		access = access + f"['{barc}']"
		
	return access
	
def construct_variable_regex(before, after, mismatches):
	"""
	Construct regex to search for variable barcodes with desired number of mismatches
	This takes the form {before}(.*){after}
	"""
	if mismatches == 0:
		return f"{before}(.*){after}"
	
	# get a regex for a mismatch in every place in before and after sequences
	befores = create_mismatches_regex([before], mismatches)
	afters = create_mismatches_regex([after], mismatches)
	
	# combine each before and after regex with (.+) in the middle
	regexes = []
	for b in befores.split("|"):
		for a in afters.split("|"):
			regexes.append(f"{b}(.*){a}")
	return "|".join(regexes)

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
		search_dict['forward_search'] = {f"{key} ({value})":value for key, value in forward_barcs.items()}
	if mismatches > 0:
		# if allowing mismatches, use regexes
		search_dict['type'] = 'constant_regex'
		search_dict['forward_search'] = {f"{key} ({value})":create_mismatches_regex([value], mismatches) for key, value in forward_barcs.items()}
		
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
		
def create_level(dict, path_list, value):
	"""
	create new level in nested dictionary, and initialise to value
	"""
	if len(path_list) == 0:
		return

	for k in path_list[:-1]:
		dict = dict[k]
	
	dict[path_list[-1]] = value
			
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
				assert len(barcode) == len(subread_forward)
				
				if barcode == subread_forward:
					matches.append(name)
					# for exact match, can only find one barcode
					break
			
			# check how many matches we found
			if not(len(matches) == 0 or len(matches) == 1):
				pdb.set_trace()
			assert (len(matches) == 0 or len(matches) == 1)
			if len(matches) == 0:
				found_barcodes.append('none')
			elif len(matches) == 1:
				found_barcodes.append(matches[0])
			
		elif set['type'] == 'constant_regex':
			matches = []
			
			# get part of read to check
			start = set['start']
			stop = set['stop']
			subread_forward = line[start:stop]
			
			# check if subread matches any of the barcodes
			for name, barcode in set['forward_search'].items():
				if re.findall(barcode, subread_forward, re.IGNORECASE):
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
			regexes =set['forward'].split("|")
			
			# match regexes
			for regex in regexes:
				match = re.findall(regex, line, re.IGNORECASE)
				for m in match:
					if m not in matches:
						matches.append(m)
						
			# check how many matches we found
			if (len(matches) == 0) or (matches[0] == "" and len(matches) == 1):
				# check if there is no insertion, or we couldn't find the regex at all
				int_type = ""
				for regex in regexes:
					# check if we found both the 'before' and 'after' sequence
					before, after = regex.split('(.*)')
					if bool(re.search(before, line, re.IGNORECASE)) and bool(re.search(after, line, re.IGNORECASE)):
						int_type = 'no_insertion'
						break
				# if we didn't find a regex above, then there wasn't a match
				if int_type == "":
					int_type = 'none'
				found_barcodes.append(int_type)

			elif len(matches) == 1:
				
				#only try to translate if there is just one match and its length is a multiple of three
				if set['trans']:
					if ( len(matches[0]) % 3 ) == 0:
						barc = str(Seq(matches[0]).translate())
					else:
						# otherwise just add brackets to indicate a nucleotide sequence
						barc = f"({matches[0]})"
				else:
					barc = matches[0]
				found_barcodes.append(barc)
				
			elif len(matches) > 1:
				found_barcodes.append('ambiguous')	

	assert len(found_barcodes) == len(search)
	return found_barcodes
	
def get_all_counts(counts, cur=()):
	"""
	Get all counts for all combinations of barcodes in the 'counts' nested dict
	https://stackoverflow.com/questions/11570499/generate-all-leaf-to-root-paths-in-a-dictionary-tree-in-python
	"""
	# if we're at the end of the tree
	if isinstance(counts, int):
		yield cur + (counts,)
	# otherwise, get all paths beyond current place in tree
	else:
		for n, s in counts.items():
			for path in get_all_counts(s, cur+(n,)):
				yield path
	
def increment_counter(counts, found_barcs):
	"""
	given a list of barcodes found in a read, and a counter dictionary
	increment the combination of found barcodes
	"""
	
	assert isinstance(counts, dict)
	assert isinstance(found_barcs, list)
	
	# first check that all levels exist in found_barcs
	for level, barc in enumerate(found_barcs):
		# check that barc exists at this level
		if barc in eval(construct_counts_access(found_barcs[:level], 'counts')):
			pass
		else:
			# need to add new level if we're not at the end of the list
			if level < len(found_barcs) - 1:
				create_level(counts, found_barcs[:level+1], {})
			# or create new entry if we are at the end of the list
			else:
				create_level(counts, found_barcs[:level+1], 0)
				
	# now increment counter
	increment_value(counts, found_barcs, 1)		
	
	
	return counts
	
def increment_value(dict, path_list, value):
	"""
	increment value in nested dictionary
	"""
	for k in path_list[:-1]:
		dict = dict[k]
	
	dict[path_list[-1]] += value
		

def interpret_bool(str):
	if str is True:
		return True
	elif str is False:
		return False
	elif str.lower() == 'true':
		return True
	elif str.lower() == 'false':
		return False
	else:
		raise ValueError(f"Couldn't understand {str} in translate field: please enter 'true' or 'false'")

def parse_barcs_yaml(args):
	"""
	parse barcodes yaml file and check that it makes sense
	"""
	
	# read barcodes yaml
	with open(args.barcodes, 'r') as stream:
		barcodes = yaml.safe_load(stream)
		
	assert isinstance(barcodes, list)
			
	print(f"found {len(barcodes)} sets of barcodes")
	names = []
	
	# check that barcodes are specified correctly
	for i in range(len(barcodes)):
		assert isinstance(barcodes[i], dict)
		assert len(barcodes[i]) == 1
		try:
			name = list(barcodes[i].keys())[0]
			names.append(name)
			# if this set is constant, check it contains a start and some sequences
			if barcodes[i][name]['type'] == 'constant':
			
				assert 'barcodes' in barcodes[i][name]
				assert 'start' in barcodes[i][name]
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
			
				assert 'before' in barcodes[i][name]
				assert 'after' in barcodes[i][name]
				
				before = barcodes[i][name]['before']
				after = barcodes[i][name]['after']
				
				# check if we shoudld translate these barcodes
				if 'translate' in barcodes[i][name]:
					translate =  interpret_bool(barcodes[i][name]['translate'])
				else:
					barcodes[i][name]['translate'] = False
					translate = False
				
				print(f"set {name} will consist of all the sequences occuring between sequences {before} and {after} in the read.  The sequences {'will' if translate else 'will not'} be translated into amino acid squences")
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
	# check for bases not recgonised
	for base in seq:
		assert base in old_chars
	return seq.translate(tab)[::-1]
		
def write_counts(outfile, counts, search):
	"""
	Write counts in recursive dictionary 'counts' as pandas data frame to file 'outfile'
	"""
	
	combinations = get_all_counts(counts)
	
	counts_df = pd.DataFrame(combinations, columns = [set['name'] for set in search] + ['count'])
	
	counts_df.to_csv(outfile, index=False, sep = '\t')
	
	

if __name__ == "__main__":
	main(sys.argv[1:])
