from flask import Flask, render_template, request, redirect, url_for, session
from flask_session import Session
from tempfile import mkdtemp
from werkzeug.utils import secure_filename
import os

from helpers import parse_form, parse_barcode_file

import pdb



UPLOAD_FOLDER = 'uploads/'
try:
	os.mkdir(UPLOAD_FOLDER)
except FileExistsError:
	pass
ALLOWED_EXTENSIONS_YAML = {'txt', 'yml', 'yaml'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


# Configure session to use filesystem (instead of signed cookies)
app.config["SESSION_FILE_DIR"] = mkdtemp()
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
Session(app)

@app.route('/')
def index():

    # Forget session
	session.clear()
	
	# prefill data if necessary
	

	return render_template('info.html')
	
def allowed_file_yaml(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS_YAML

@app.route('/', methods=['POST'])
def barcodes():

	# check if there was a file uploaded
	file = request.files['file']

	# if no file
	if file.filename != "":
		
		# check file is allowed type
		if not allowed_file_yaml(file.filename):
			render_template('info.html', error = "Incorrect file type")			
		
		# save yaml file
		filepath = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(file.filename))
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
	
	print(session)
	
	if 'barcodes' not in session:
		return redirect(url_for('index'))
	
	
	

	return render_template("files.html")
    
if __name__ == "__main__":
	app()