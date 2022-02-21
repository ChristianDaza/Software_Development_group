# Simple example of a web application
# pp05: Now with a form to enter protein name
from flask import Flask, render_template, url_for, redirect, abort, send_from_directory

######### Chris' imports #######
import pandas as pd
import os
from flask_apscheduler import APScheduler
import atexit
import time
################################################

# import libraries needed create and process forms
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import InputRequired



# create a flask application object
app = Flask(__name__)


######### Needed for  Chris'  fucntions #######
# Set configuration values
class Config:
    SCHEDULER_API_ENABLED = True

# Initialize scheduler
scheduler = APScheduler()

################################################



# we need to set a secret key attribute for secure forms
app.config['SECRET_KEY'] = 'change thisunsecure key'

######### Needed for  Chris'  fucntions #######
# Path to download file
app.config['path_stat_results_folder'] = os.path.join(os.getcwd(), 'dynamic')
################################################

# tell code where to find protein information
protein_table_filename = 'Book1.csv'


# create a class to define the form
class QueryForm(FlaskForm):
    protein_name = StringField('Enter a valid UniProt protein name:', validators=[InputRequired()])
    submit = SubmitField('Submit')


# define the action for the top level route
@app.route('/', methods=['GET', 'POST'])
def index():
    # this route has been updated to use a template containing a form
    form = QueryForm()  # create form to pass to template
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        #save_input =  request.form.get(protein_name)
        print('\n\n\n' + protein_name + '\n\n\n')
        return redirect(url_for('protein', protein_name=protein_name))
    return render_template('index_page.html', form=form, protein_name=protein_name)


# define a route called 'protein' which accepts a protein name parameter
@app.route('/protein/<protein_name>')
def protein(protein_name):
    print(protein_name)
    # load protein data from TSV file into pandas dataframe with protein name as index
    print("This part works")
    df = pd.read_csv('Book1.csv')
    # ensure name is in capital letters
    # try to extract row for specified protein
    row = df.query("ID == '%s'" % protein_name)
    ########## Chri's modifications to search fucntion ###############

    # Check that search results are not empty and that the folder where the results are going to be save exist, if not create it
    if row.empty == False and not os.path.exists("dynamic"):
        os.makedirs("dynamic")
        row.to_csv(path_or_buf=os.path.join(os.getcwd(), 'dynamic', "{}.csv".format(protein_name)))
        return render_template('protein_view.html', tables=[row.to_html(classes='data')], titles=row.columns.values)

    elif row.empty == False:
        row.to_csv(path_or_buf=os.path.join(os.getcwd(), 'dynamic', "{}.csv".format(protein_name)))
        return render_template('protein_view.html', tables=[row.to_html(classes='data')], titles=row.columns.values)

    else:
        abort(404, "No data found, please try with a different name")

################################ Chris' functions  ######################################

#  Look into the dynamic folder and create a list of save files

@app.route('/display_downloads/')
def display_downloads():
    """ Takes the names of the statistics results  saved as csv and appends it to a list for HTML display"""
    results_stats_filenames = os.listdir(app.config['path_stat_results_folder'])
    results_stats_filenames_list = []
    # Append each statistics results file into the list
    try:
        for file in results_stats_filenames:
            results_stats_filenames_list.append(file)
        return render_template('file_downloads.html',  results_stats_filenames_list = results_stats_filenames_list)
    except:
        abort(404, "There was a problem displaying your results, please try again")




@app.route('/get_csv/<filename>', methods=["GET"])
def get_csv(filename):
    """ Download statistics results csv file/s into users computer and  delete them after"""
    # Download and then delete  the  user's chosen saved statistics results if the file is not empty
    try:
        results_stats_filenames = os.stat(os.path.join(os.getcwd(), 'dynamic', "{}".format(filename))).st_size
        if results_stats_filenames > 0:
            return send_from_directory(app.config['path_stat_results_folder'], filename), \
                   os.remove(os.path.join(os.getcwd(), 'dynamic', "{}".format(filename))),
    except:
        abort(404, "Empty file, please select another file")



@scheduler.task('interval', id='del_files', seconds=60, misfire_grace_time=900)
def del_files():
    """ Delete the saved statistics  results  csv files after 30 minutes """
    current_time = time.time()
    size_stats_results_folder = os.stat(app.config['path_stat_results_folder']).st_size

    # if folder of save statistics results is not empty and the csv files inside are less tha a day old, delete them
    if size_stats_results_folder > 0:
        for stat_results_filename in os.listdir(app.config['path_stat_results_folder']):
            if os.path.getmtime(os.path.join(app.config['path_stat_results_folder'], stat_results_filename)) < current_time   -1:
                if os.path.isfile(os.path.join(app.config['path_stat_results_folder'], stat_results_filename)):
                    print(stat_results_filename)
                    os.remove(os.path.join(app.config['path_stat_results_folder'], stat_results_filename))
    else:
        pass






# Deletion fucntion




# start the web server
if __name__ == '__main__':
    ######### Needed for  Chris'  fucntions #######
    scheduler.start()
    atexit.register(lambda: scheduler.shutdown())
    ################################################
    app.run(debug=True)

