from flask import render_template,abort,Blueprint, Flask
from jinja2 import TemplateNotFound

app=Flask(__name__)

@app.route("/")
def showCrisprCBEI():
    try:
        return render_template('index.html')
    except TemplateNotFound:
        abort(404)