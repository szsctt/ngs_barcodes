import re
import yaml
import os
import pdb

from werkzeug.utils import secure_filename

from setup import ALLOWED_EXTENSIONS_YAML


def parse_barcode_file(filename):
	"""
	Parse barcode file to generate representation of the barcode information
	Check that this looks to be the correct format
	"""
	
	with open(filename, "r") as stream:
		try:
			barcodes = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			return None, "Not valid YAML"
			
	# barcodes must be a list
	if not isinstance(barcodes, list):
		return None, "Can't identify sets of barcodes in uploaded file"
	
	if len(barcodes) == 0:
		return None, "Can't identify any sets of barcodes uploaded file"		
	
	# check each barcode set
	for set in barcodes:
		
		# if set isn't a dictionary
		if not isinstance(set, dict):
			return None, f"Set number {i+1} must contain key-value pairs"
		
		# must have only one key/value pair, otherwise it's not clear what to do with the others
		if len(set) != 1:
			return None, f"Set number {i+1} must only have one entry"		
		
		set_name = list(set.keys())[0]
		print(f"checking '{set_name}'")			
		
		# check set values
		set_info = set[set_name]
		
		# must have a key 'type' which has value 'constant' or 'variable'
		if not 'type' in set_info:
			return None, f"Set '{set_name}' must have a key 'type' with value 'constant' or 'variable'"			
			
		# check if set conforms to specifications for set type
		if set_info['type'] == 'constant':
			err = check_constant_set(set_info, set_name)
		elif set_info['type'] == 'variable':
			err = check_variable_set(set_info, set_name)
		else:
			err = f"Set '{set_name}' must have a key 'type' with value 'constant' or 'variable'"
			
		if err:
			return None, err
	
	return barcodes, None

def check_constant_set(set, set_name):
	"""
	Check if a set of barcodes is valid
	Return an error message if not valid, otherwise return None if set is valid
	"""
	
	# must have key 'start' with a  positive integer value
	if 'start' not in set.keys():
		return f"Constant set '{set_name}' must have a key 'start' which specifies where the barcode starts in the read"
	
	if not isinstance(set['start'], int):
		return f"Key 'Start' in constant set '{set_name}' must have a positive interger value"
	
	if set['start'] < 0:
		return f"'Start' key for constant set '{set_name}' must have a positive interger value"
		
	# must have a key 'barcodes' containing a dict of valid barcodes, all with the same length and composed of a/c/g/t
	if 'barcodes' not in set.keys():
		return f"Constant set '{set_name}' must have a key 'barcodes' with contains a dictionary of the barcodes to look for"
		
	if not isinstance(set['barcodes'], dict):
		return f"Key 'barcodes' in constant set '{set_name}' must contains a dictionary of the barcodes to look for"
		
	if len(set['barcodes']) == 0:
		return f"Key 'barcodes' in constant set '{set_name}' must contains a dictionary of the barcodes to look for"
		
	# check each barcode
	seq_len = len(list(set['barcodes'].values())[0])
	for name, seq in set['barcodes'].items():
		
		# check sequence composed of a/c/g/t
		if not check_seq(seq):
			return f"Barcode '{name}' in set '{set_name}' must be composed of only 'A', 'C', 'G' and 'T'"
		
		# check same length as first sequence
		if not len(seq) == seq_len:
			return f"Barcodes in set '{set_name}' must all be the same length"
	
	return None
	
def check_variable_set(set, set_name):
	"""
	Check if a variable set of barcodes is valid
	Return an error message if not valid, otherwise return None if set is valid
	"""
	
	# must have key 'before' composed of a/c/g/t
	if 'before' not in set.keys():
		return f"Constant set '{set_name}' must have a key 'before' which specifies the sequence preceeding the variable part of the read"
	
	if not check_seq(set['before']):
		return f"Key 'before' in constant set '{set_name}' must be a nucleotide sequence (composed of only 'A', 'C', 'G' and 'T')"
	
	# must have key 'after' composed of a/c/g/t
	if 'after' not in set.keys():
		return f"Constant set '{set_name}' must have a key 'after' which specifies the sequence preceeding the variable part of the read"
	
	if not check_seq(set['after']):
		return f"Key 'after' in constant set '{set_name}' must be a nucleotide sequence (composed of only 'A', 'C', 'G' and 'T')"	
	
	return None

def check_seq(seq):
	"""
	Check that a string looks like a DNA sequence (composed of only A/C/G/T)
	"""
	
	# sequence must be nonzero length
	if len(seq) == 0:
		return False
	
	return all([let.lower() in {'a', 'c', 'g', 't'} for let in seq])	

def parse_form(form):
	"""
	Form returns flat data format - parse form input to reproduce barcodes config format
	"""
	
	print(form)
	
	set_names = [key for key in form.keys() if re.search("set_[a-zA-Z]+_name", key)]
	set_letters = [re.search("set_([a-zA-Z]+)_name", key).group(1) for key in set_names]
	
	yaml = []
	errors = []
	
	for let in set_letters:
	
		# collect all key/value pairs from this set
		set = {key:form[key] for key in form.keys() if re.search(let, key)}
		
		# get set name
		name = form[f"set_{let}_name"]
	
		# is set constant or variable?
		if f"set_{let}_before" in set:
		
			# parse set
			yaml_set = parse_variable(set, let)
			
			# check set is valid
			err = check_variable_set(yaml_set, name)

		else:
		
			# parse set
			yaml_set = parse_constant(set, let)
			
			# check set is valid
			err = check_constant_set(yaml_set, name)
			
		if err:
			errors.append(err)
		
		yaml.append({name: yaml_set})
		
	
	return yaml, errors

def num_to_str(num):
	"""
	Convert a number to a string, where 0 is 'a', 25 is 'z', 26 is 'za' and 27 is 'zb'
	"""
	
	num_zs = num // 26
	last_let = num % 26
	
	return 'z' * num_zs + chr(97 + last_let)
	
def parse_constant(set, set_letter):
	"""
	Convert flat dict to correct format for a constant set of barcodes
	"""
	
	yaml_set = {'type': 'constant'}
	
	# get set name
	name_key = [key for key in set if re.search(f"set_{set_letter}_name", key)][0]
	yaml_set['name'] = set[name_key]
	
	# get set position
	pos_key = [key for key in set if re.search(f"set_{set_letter}_pos", key)][0]
	yaml_set['start'] = set[pos_key]
	
	# get set barcodes
	yaml_set['barcodes'] = {}
	i = 0
	while f"{set_letter}_{num_to_str(i)}_name" in set:
		yaml_set['barcodes'][set[f"{set_letter}_{num_to_str(i)}_name"]] = set[f"{set_letter}_{num_to_str(i)}_seq"]
		i += 1
		
	
	return yaml_set
		
def parse_variable(set, set_letter):
	"""
	Convert flat dict to correct format for a variable set of barcodes
	"""
	
	yaml_set = {'type': 'variable'}

	# get set before seq
	yaml_set['before'] = set[f"set_{set_letter}_before"]
	
	# get set after seq
	yaml_set['start'] = set[f"set_{set_letter}_after"]
	
	# get if we should translate this set
	yaml_set['translate'] = set[f"{set_letter}_translate"]
	
	return yaml_set
	
def allowed_file_yaml(filename):
	"""
	Check if filename extension is allowed for YAML
	"""
	return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS_YAML

def save_read_files(request_files, session):
	"""
	Save R1 and R2 files input by user
	"""
	
	R1_file_keys = [key for key in request_files if re.search("R1_\d+", key)]
	session['fastq_files'] = []
	
	
	for i, R1 in enumerate(R1_file_keys):
		
		# check for no file uploaded
		if request_files[R1].filename == "":
			del session['fastq_files']
			return "Please specify all files"
			
		session['fastq_files'].append([])
		
		# save R1
		filepath = os.path.join(session['folder'], secure_filename(request_files[R1].filename))
		request_files[R1].save(filepath)
		session['fastq_files'][i].append(filepath)
			
		# save R2
		R2 = re.sub("R1_", "R2_", R1)	
				
		# check for no file uploaded
		if request_files[R2].filename == "":
			del session['fastq_files']
			return "Please specify all files"
			
		
		filepath = os.path.join(session['folder'], secure_filename(request_files[R2].filename))
		request_files[R2].save(filepath)
		session['fastq_files'][i].append(filepath)
		


def create_filesets(request_form, session):
	"""
	Save information about which barcodes go with which files to session
	"""
	
	session['fastq_barcode_assoc'] = []
	
	fastqset_names = [key for key in request_form.keys() if re.search("fastqset_name_[0-9]+", key)]
	fastqset_nums = [int(re.search("fastqset_name_([0-9]+)", key).group(1)) for key in fastqset_names]
	
	for num in fastqset_nums:
	
		session['fastq_barcode_assoc'].append({})
		
		# get info from form for this dataset
		name =  request_form[f"fastqset_name_{num}"]
		f1 = request_form[f"R1_{num}"]
		f2 = request_form[f"R1_{num}"]
		
		if f1 == "" or f2 == "":
			return "Please specify all read files"
	
		# if a common part of the filenames was specified
		if name != "":
			suffix1, suffix2 = check_common_part(name)
			
			if suffix1 == "" or suffix2 == "":
				return f"Couldn't find name {name} in read file names {f1} or {f2}"
			
		else:
			# get common part of the name and suffixes
			name, suffix1, suffix2 = get_common_part(f1, f2)
			
		
			if sample == "":
				return f"Read file names {f1} and {f2} do not contain a common part. Are you sure these are correctly paired?"
		
		# get adapters
		
		
		# get lengths
		

		
		# add info to session
		session['fastq_barcode_assoc'][num]['name'] = name
		session['fastq_barcode_assoc'][num]['R1_suffix'] = suffix1
		session['fastq_barcode_assoc'][num]['R2_suffix'] = suffix2

	
		pdb.set_trace()
		print(nums)

def check_common_part(part, f1, f2):
	"""
	Check that the common part of a filename occurs at the start of f1 and f2.  Return a tuple of the f1 and f2 suffixes (non-common part of the filenames)
	"""
	
	regex = f"^{re.escape(part)}(.+)$"
	
	m1 = re.match(regex, f1)
	if not m1:
		return "", ""
	
	m2 = re.match(regex, f2)
	if not m2:
		return "", ""
		
	return m1.group(1), m2.group(2)
	

def get_common_part(f1, f2):
	"""
	Get the common part at the start of filenames f1 and f2
	"""
	common = ""
	for i, let in enumerate(f1):
		if f2[i] == let:
			common = common + let
		else:
			break
			
	if i >= len(f1) or i >= len(f2):
		return "", "", ""
			
	suffix1 = f1[i:]
	suffix2 = f2[i:]
		
			
	return common, suffix1, suffix2
	
	