# python imports
import os
import shutil
import tempfile
import secrets

# flask imports
from flask import Flask, render_template, request, redirect, url_for, session, flash, send_file
#from flask_session import Session
from werkzeug.utils import secure_filename

# redis imports
import redis
from rq import Queue
from rq.job import Job

# helpers
from helpers import parse_form, parse_barcode_file, allowed_file_yaml, save_read_files, create_filesets, run_snakemake
from setup import UPLOAD_FOLDER

import pdb

# make folder for uploads, if it doesn't already exist
try:
	os.mkdir(UPLOAD_FOLDER)
except FileExistsError:
	pass

# initialize flask app
app = Flask(__name__)

# configure upload folder
app.config['UPLOAD_FOLDER'] = tempfile.mkdtemp()


# secret key for sessions
app.secret_key = secrets.token_hex()

print(f"using folder {app.config['UPLOAD_FOLDER']} for uploads")

# redis initilzation
redis_conn = redis.Redis(host='redis', port=6379)
q = Queue(connection=redis_conn, default_timeout=7200)

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

	# make sure we've already got some barcodes and files to work with
	check = check_session(barcodes = True, files = True)
	if check is not None:
		return check

	# check for a redis job
	if 'job' not in session:
	
		# enqueue job
		# can't pickle session, so make copy of items
		job = q.enqueue(run_snakemake, dict(session.items()))
		
		session['job'] = job.id
		
		return render_template("results.html", finished=False)
	
	else:
	
		# get job status
		job = Job.fetch(session['job'], connection = redis_conn)
		status = job.get_status()
		
		print(f'Job status: {status}')
		
		# if something went wrong
		if status in ("canceled", "failed"):

			flash('Something went wrong!  Please try again')
			return redirect(url_for('index'))
			
		# if not finished yet
		if status in ("queued", "started", "deferred", "scheduled"):
		
			return render_template("results.html", finished=False)
			
		else:


			#send_from_directory(job.result)
			return render_template("results.html", finished=True)
		
@app.route('/return-files/')
def return_files_tut():
	
	job = Job.fetch(session['job'], connection = redis_conn)

	try:
		return send_file(job.result)
	except Exception as e:
		return str(e)	

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