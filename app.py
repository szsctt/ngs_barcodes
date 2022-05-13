import os
import shutil
import tempfile
import secrets

from flask import Flask, render_template, request, redirect, url_for, session
from flask_session import Session
from werkzeug.utils import secure_filename

from helpers import parse_form, parse_barcode_file, allowed_file_yaml, save_read_files, create_filesets, run_snakemake
from setup import UPLOAD_FOLDER

import pdb


try:
	os.mkdir(UPLOAD_FOLDER)
except FileExistsError:
	pass


app = Flask(__name__)

# configure upload folder
app.config['UPLOAD_FOLDER'] = tempfile.mkdtemp()


# secret key for sessions
app.secret_key = secrets.token_hex()

print(f"using folder {app.config['UPLOAD_FOLDER']} for uploads")

@app.route('/')
def index():

	print("index, get")
	# clear temporary directory, if there is one
	if 'folder' in session:
		shutil.rmtree(session['folder'])

    # Forget session
	session.clear()
	
	# make a new temporary folder for this session
	folder = tempfile.mkdtemp(dir = app.config['UPLOAD_FOLDER'])
	session['folder'] = folder
	
	return render_template('info.html')

@app.route('/', methods=['POST'])
def barcodes():

	print("index, post")

	# check if there was a file uploaded
	file = request.files['file']

	# if no file
	if file.filename != "":
		
		# check file is allowed type
		if not allowed_file_yaml(file.filename):
			render_template('info.html', error = ["Incorrect file type"])			
		
		# save yaml file
		filepath = os.path.join(session['folder'], secure_filename(file.filename))
		file.save(filepath)
		
		# read barcodes from file
		barcodes, err = parse_barcode_file(filepath)

		if err:
			err = [err]
		else:
			err = []
    
    # otherwise, check for entered barcodes in form
	else:
		barcodes, err = parse_form(request.form)

	# return page with error, if there was one
	if len(err) > 0:
		return render_template('info.html', errors = err)
    
	# save barcodes to session
	session['barcodes'] = barcodes

	return redirect(url_for('files'))
	
@app.route("/files")
def files():

	print("files, get")
	
	check = check_session(barcodes=True)
	if check is not None:
		return check
		
	barcode_names = [list(set.keys())[0] for set in session['barcodes']]

	return render_template("files.html", barcode_names = barcode_names)
	
@app.route('/files', methods=['POST'])
def submit_files():
	
	# save read files
	err = save_read_files(request.files, session)
	barcode_names = [list(set.keys())[0] for set in session['barcodes']]
	if err:
		return render_template("files.html", error = err, barcode_names=barcode_names)
		
	# check rest of form info and save to session
	err = create_filesets(request.form, session)
	if err:
		return render_template("files.html", error = err, barcode_names=barcode_names)	
	
	# redirect to results
	return redirect(url_for('results'))


@app.route('/results')
def results():


	pdb.set_trace()
	check = check_session(barcodes = True, files = True)
	if check is not None:
		return check
		
	# if we haven't already started the analysis
	if 'proc' not in session:
		run_snakemake(session)
		

	if session['proc'].poll() is None:
		return render_template("results.html")
		
	

def check_session(barcodes = False, files = False):
	
	if barcodes:
	
		if 'barcodes' not in session:
			return redirect(url_for('index'))
		
		if len(session['barcodes']) == 0:
			return redirect(url_for('index'))	
	
	if files:
		
		if 'fastq_files' not in session:
			return redirect(url_for('files'))
		

    
if __name__ == "__main__":
	app()