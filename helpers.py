import re
import pdb

def parse_form(form):
	"""
	Parse form input to re-produce barcodes config format
	"""
	
	print(form)
	
	set_names = [key for key in form.keys() if re.search("set_[a-zA-Z]+_name", key)]
	set_letters = [re.search("set_([a-zA-Z]+)_name", key).group(1) for key in set_names]
	
	yaml = []
	
	for let in set_letters:
	
		# collect all key/value pairs from this set
		set = [key for key in form.keys() if re.search(let, key)]
	
		# is set constant or variable?
		if f"set_{let}_before" in set:
			yaml.append(parse_variable(set))
		else:
			yaml.append(parse_constant(set))
		
	
#	pdb.set_trace()
	
def parse_constant(set):
	"""
	Convert flat dict to correct format for a constant set of barcodes
	"""
	
	yaml_set = {'type': 'constant'}
	
	# get set name
	name_key = [key for key in set if re.search("set_[a-zA-Z]+_name", key)][0]
	yaml_set['name']
	
	# get 
	
	
	return yaml_set
	
	
def parse_variable(set):
	"""
	Convert flat dict to correct format for a variable set of barcodes
	"""
	
	yaml_set = {'type': 'variable'}
	
	
	return yaml_set