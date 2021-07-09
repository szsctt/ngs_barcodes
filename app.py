from flask import Flask, render_template, request, redirect, url_for, session
from flask_session import Session
from tempfile import mkdtemp

from helpers import parse_form

import pdb

app = Flask(__name__)


# Configure session to use filesystem (instead of signed cookies)
app.config["SESSION_FILE_DIR"] = mkdtemp()
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
Session(app)

@app.route('/')
def index():

    # Forget session
	session.clear()

	return render_template('info.html')

@app.route('/', methods=['POST'])
def barcodes():
# 
# 	# make sure all files were defined
# 	if all([uploaded_file.filename != '' for uploaded_file in  request.files.getlist('file')]):
# 		for uploaded_file in request.files.getlist('file':
# 			uploaded_file.save(uploaded_file.filename)
# 

    # Forget session
	session.clear()
    
	# save barcodes to session for now
	session['barcodes'] = parse_form(request.form)

	return redirect(url_for('files'))
	
@app.route("/files")
def files():
	
	print(session)
	
	if 'barcodes' not in session:
		return redirect(url_for('index'))
	
	
	

	return render_template("files.html")
    
if __name__ == "__main__":
	app()