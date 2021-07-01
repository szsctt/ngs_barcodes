from flask import Flask, render_template, request, redirect, url_for
import pdb

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('info.html')

@app.route('/', methods=['POST'])
def upload_file():

	print(request.form)
# 
# 	# make sure all files were defined
# 	if all([uploaded_file.filename != '' for uploaded_file in  request.files.getlist('file')]):
# 		for uploaded_file in request.files.getlist('file':
# 			uploaded_file.save(uploaded_file.filename)
# 
	return redirect(url_for('index'))
    
if __name__ == "__main__":
	app()