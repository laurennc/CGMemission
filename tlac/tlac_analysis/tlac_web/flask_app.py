"""Run tlac_web embedded in a flask app.


How to run it:
cd to_a_directotry_where_there_are_tlac_files
bokeh-server --script PATH_TO_HERE/tlac_web/tlac_web.py
(new terminal)
cd to_a_directotry_where_there_are_tlac_files
python PATH_TO_HERE/tlac_web/flask_app.py

For information on usage from a user perspective see the README.md in this
directory.


Todo:
- Automatic reloading of uploaded file list
- Change order of line info and db info
- Chisq

"""
import logging
LOG_FORMAT = '[flask_app] %(asctime)-15s %(message)s'
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

from bokeh.pluginutils import app_document
import markdown

from flask import Flask, render_template, url_for, send_file, request,\
    redirect, Markup
from werkzeug import secure_filename

from tlac_web import TlacWebApp

import os, sys, glob, time
import tarfile

import numpy as np
import uuid

## Run time arguments
public_mode = False
data_dir = None
upload_dir = None
bokeh_url = "http://localhost:5006"
port = 5050

for i,carg in enumerate(sys.argv):
    if carg == "--public":
        print "Running in public mode"
        public_mode = True
    if carg == "--data_dir":
        print "Enabling data download."
        data_dir = os.path.realpath(sys.argv[i+1]) + "/"
        print "Saving temporary data in " + data_dir
    if carg == "--upload_dir":
        print "Enabling data upload."
        upload_dir = os.path.realpath(sys.argv[i+1]) + "/"
        print "Saving temporary data in " + upload_dir
    if carg == '--bokeh_url':
        bokeh_url = sys.argv[i+1]
    if carg == '--port':
        port = int(sys.argv[i+1])
        
print "Looking for bokeh-server on " + bokeh_url

app = Flask('tlac_web')

ALLOWED_EXTENSIONS = set(['txt', 'dat','tsv'])

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 0.1 * 1024 * 1024 # 0.1MB
app.config['UPLOAD_FOLDER'] = upload_dir
app.config['PORT'] = port

################################################################################

def get_data_prefix(sid, mydir):        
    if mydir is not None:
        if "*" in sid:
            print "Problem encountered. sid=",sid
            sid = make_session_id()
        data_prefix = mydir + sid
    else:
        data_prefix = None

    return data_prefix

def make_session_id():
    return time.strftime("%Y-%m-%d_%H%M_") + uuid.uuid4().hex

    
@app_document("tlac_web", bokeh_url)
def make_tlac_web_applet(session_id):
    app = TlacWebApp.create(data_prefix = get_data_prefix(session_id, data_dir),
                            upload_prefix = get_data_prefix(session_id,
                                                            upload_dir))
    return app

@app.route("/tlac_web")
def applet():
    session_id = make_session_id()
    logging.info("New tlac_web session " + session_id)
    applet = make_tlac_web_applet(session_id)
    if data_dir is None:
        dl_link = ""
    else:
        dl_link = url_for('download_data', session_id = session_id)

    if upload_dir is None:
        ul_link = ""
    else:
        ul_link = url_for('upload_data', session_id = session_id)
    
    return render_template(
        "tlac_web_01.html",
        app_url = bokeh_url + "/bokeh/jsgenerate/VBox/TlacWebApp/TlacWebApp",
        app_tag = applet._tag,
        dl_link = dl_link,
        ul_link = ul_link
    )

@app.route("/tlac_web/download_data/<session_id>")
def download_data(session_id):
    logging.info("Session " + session_id + " downloads data.")

    if session_id is None or session_id == "":
        return ""
    data_prefix = get_data_prefix(session_id, data_dir)
    filelst = glob.glob(data_prefix + "*.dat")

    if len(filelst) == 0:
        return "No files generated."\
            "Please use the back button of your browser and generate data "\
            "files by pressing the \"Save\" button in tlac_web."

    tarfn = data_prefix + ".tar.gz"
    tar = tarfile.open(tarfn, "w:gz")
    for name in filelst:
        tar.add(name, arcname=os.path.basename(name))
    tar.close()

        
    return send_file(tarfn, as_attachment=True)

@app.route("/tlac_web/upload_data/<session_id>", methods=['GET', 'POST'])
def upload_data(session_id):
    logging.info("Session " + session_id + " uploads data.")
    if session_id == "":
        return ""

    warn_msg = ""
    if request.method == 'POST':
        file = request.files['file']
        if file:
            if allowed_file(file.filename):
                #filename = secure_filename(file.filename)
                filename = secure_filename(session_id + file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                return redirect(url_for('uploaded_file',
                                        filename=filename))
            else:
                warn_msg = "Please select data with valid extension."
                warn_msg += "(" + ", ".join(ALLOWED_EXTENSIONS) + ")"
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <form action="" method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    <strong>''' + warn_msg + '''</strong>'''

    
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[-1] in ALLOWED_EXTENSIONS


@app.route('/tlac_web/uploads/<filename>')
def uploaded_file(filename):
    try:
        dat = np.loadtxt(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    except:
        msg = "Error: Please select whitespace-seperated "\
              "data file."
        os.remove(filename)
        return msg

    shape = np.shape(dat)
    if 2 in shape:        
        msg = "Upload successful!<br />Uploaded file with (nrows, ncols)="+\
              str(shape) + "<br />Change any slider to update the file list."
    else:
        msg = "Please upload file with two rows or columns."
        os.remove(filename)

    return msg

@app.route('/tlac_web/doc',methods=['GET'])
def doc_as_html():
    """Parses markdown of README.md in current folder and outputs it.

    Info:  http://flask.pocoo.org/snippets/19/
    """
    logging.info("Docs are opened.")
    
    fn = os.path.dirname(os.path.realpath(__file__)) + "/README.md"
    
    if not os.path.isfile(fn):
        return "(no README found)"
    f = open(fn,"r")
    doc = "".join(f.readlines())
    f.close()

    return Markup(markdown.markdown(doc))
    
    

if __name__ == "__main__":
    print "Running tlac_web embedded in flask at port " +str(app.config['PORT'])
    if not public_mode:
        app.debug = True
        app.run(port=app.config['PORT'])
    else:
        app.run(host='0.0.0.0', port=app.config['PORT'])

