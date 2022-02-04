from flask import Flask, render_template
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import  DataRequired

app = Flask(__name__)
#app.config['SECRET_KEY'] = "super secret key"

#create a form class

#class rsID_Form(FlaskForm):
#    rsID = StringField("Enter rsID value", validators=[DataRequired()])
#    Submit = SubmitField("Submit")

@app.route('/')
def home():
    return render_template('base.html')

@app.route('/search')
def search():
    return render_template('home.html')

@app.route('/About')
def About():
    return render_template('About.html')

@app.route('/Tutorial')
def Tutorial():
    return render_template('Tutorial.html')


if __name__ == '__main__':
    app.run(debug=True)
