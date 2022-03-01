
from flask import Flask, render_template, url_for, redirect, request, flash, send_file
import pandas as pd
import os
from flask_apscheduler import APScheduler
import atexit
import time

# import libraries needed create and process forms
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, IntegerField
from wtforms.validators import InputRequired
from flask_bootstrap import Bootstrap
from stats import nucleotide_diversity, homozygosity, tajima_D, fst, nucleotide_diversity_plot, homozygosity_plot, tajima_D_plot, fst_plot

# create a flask application object
app = Flask(__name__)
# initialise flask_bootstrap
Bootstrap(app)

# we need to set a secret key attribute for secure forms
app.config['SECRET_KEY'] = 'change thisunsecure key'

# Set the path for the  saved  statistical results file
app.config['path_Results.txt'] = os.path.join(os.getcwd(), 'Downloads', 'Results.txt')

# Set the path to the download folder
app.config['path_Download_folder'] = os.path.join(os.getcwd(), 'Downloads')

# Set configuration for app scheduler
class Config:
    SCHEDULER_API_ENABLED = True

# Initialize scheduler
scheduler = APScheduler()

# create classes to define the forms
class rsIDForm(FlaskForm):
    #string search form with a validator (something has to be input)
	rsID = StringField('rsID', validators=[InputRequired()])
	submit = SubmitField(label='Submit')


class posForm(FlaskForm):
    # two integer search forms both with validators
	position1 = IntegerField('Enter the first position:', validators=[InputRequired()])
	position2 = IntegerField('Enter the second position:', validators=[InputRequired()])
	submit2 = SubmitField('Submit')

class geneForm(FlaskForm):
    gene = StringField('Enter a valid HGNC accepted Gene name or alias:', validators=[InputRequired()])
    submit3 = SubmitField('Submit')


# class to store data in global environment to allow for access of information between multiple routes
class DataStore():
    pos1=None
    pos2=None
    genepos1=None
    genepos2=None
    gname=None
data = DataStore()

# read database
dff = pd.read_csv('database.csv.gz', compression="gzip")

# define the top level route
@app.route('/')
def home():
    # returns template specified
    return render_template('home-page.html')

@app.route('/about')
def about():
    return render_template('about.html')

# defining the action of the 'rsID' route
@app.route('/rsID', methods=['GET','POST'])
def rsID():
    form = rsIDForm()  # create form to pass to template
    rs_ID = None
    if form.validate_on_submit():
        rs_ID = form.rsID.data  # assign rs_ID as the user input
        return redirect(url_for('rsID_results', rsID = rs_ID))   # redirect to rsID results function with rsID as a parameter
    return render_template('rsID.html', form=form, rsID= rs_ID) # return to rsID page

# defining the action of the 'position' route
@app.route('/position', methods=['GET','POST'])
def pos():
    form = posForm()  # create form to pass to template
    position_1 = None
    position_2 = None
    if form.validate_on_submit():
        position_1 = form.position1.data
        position_2 = form.position2.data
           # redirect to position results function with position1 and position2 as parameters
        return redirect(url_for('pos_results' ,position1 = position_1, position2=position_2))
    flash('Please enter a valid numeric range.')    # if no valid input message is flased and template rendered
    return render_template('position.html', form=form, position1 = position_1, position2=position_2)



# defining the action of the 'Gene_name' route
@app.route('/Gene_name', methods=['GET','POST'])
def gene():
    form = geneForm()  # create form to pass to template
    gene = None
    value = request.form.get('check')
    if form.validate_on_submit():
        gene = form.gene.data
        # redirect to gene results function with gene as a parameter
        return redirect(url_for('gene_results', gene = gene))
    return render_template('gene.html', form=form, gene = gene)




# define rsIDresults route that takes an rsID parameter
@app.route('/rsIDresults/<rsID>')
def rsID_results(rsID):
    # query database for an rsID using pandas.query and store results in a new dataframe
    df3 = dff.query("ID == '%s'" % rsID)
    if df3.empty == False:
        # manipulate dataframe into multiple smaller dfs for each population
        X = df3.iloc[:, 0:7]
        GIH = df3.iloc[:, [1, 2, 7, 8, 9, 10, 11]]
        GWD = df3.iloc[:, [1, 2, 12, 13, 14, 15, 16]]
        JPT = df3.iloc[:, [1, 2, 17, 18, 19, 20, 21]]
        PUR = df3.iloc[:, [1, 2, 22, 23, 24, 25, 26]]
        TSI = df3.iloc[:, [1, 2, 27, 28, 29, 30, 31]]
        # hyperlink the position of each rsID to their respective regions in ensembl
        X['POS'] = X['POS'].apply(lambda x: f'<a href="https://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;r=22:{x}">{x}</a>')
        # turning dataframes of each population into html tables using df.to_html
        return render_template('rsID_results.html', tables=[df3.to_html(classes='data')], titles=df3.columns.values, rsID=rsID, table1=[X.to_html(classes='table text-center bottommarg', index=False, escape=False)], title1=X.columns.values,
                           table2=[GIH.to_html(classes='table text-center bottommarg', index=False, escape=False)], title2=GIH.columns.values,
                           table3=[GWD.to_html(classes='table text-center bottommarg', index=False, escape=False)], title3=GWD.columns.values,
                           table4=[JPT.to_html(classes='table text-center bottommarg', index=False, escape=False)], title4=JPT.columns.values,
                           table5=[PUR.to_html(classes='table text-center bottommarg', index=False, escape=False)], title5=PUR.columns.values,
                           table6=[TSI.to_html(classes='table text-center bottommarg', index=False, escape=False)], title6=TSI.columns.values)
    else:
        # flash message if no results in dataframe rather than showing an empty one
        flash('No results found.')
        # return to route currently on (rsID)
        return redirect('/rsID')



@app.route('/pos_results/<position1>/<position2>', methods=['GET','POST'])
def pos_results(position1, position2):
    # query database for rows where POS is higher that position 1 creating the rows df
    rows = dff.query("POS > %s" % position1)
    # query rows dataframe for POS less than position 2
    rows2 = rows.query("POS < %s" % position2)
    if rows2.empty==False:
        data.pos1 = int(position1)
        data.pos2 = int(position2)
        # determine how many results we obtain
        hits = rows2.shape[0]
        # hyperlink the ID of each row to their respective rsIDresults page
        rows2['ID'] = rows2['ID'].apply(lambda x: f'<a href="/rsIDresults/{x}">{x}</a>')
        # manipulate dataframe into multiple smaller df's for each population
        X = rows2.iloc[:, 0:7]
        GIH = rows2.iloc[:, [1, 2, 7, 8, 9, 10, 11]]
        GWD = rows2.iloc[:, [1, 2, 12, 13, 14, 15, 16]]
        JPT = rows2.iloc[:, [1, 2, 17, 18, 19, 20, 21]]
        PUR = rows2.iloc[:, [1, 2, 22, 23, 24, 25, 26]]
        TSI = rows2.iloc[:, [1, 2, 27, 28, 29, 30, 31]]
        # turning dataframes of each population into html tables using df.to_html
        return render_template('position_results.html', hits=hits, table1=[X.to_html(classes='table text-center bottommarg', index=False, escape=False)], title1=X.columns.values,
                           table2=[GIH.to_html(classes='table text-center bottommarg', index=False, escape=False)], title2=GIH.columns.values,
                           table3=[GWD.to_html(classes='table text-center bottommarg', index=False, escape=False)], title3=GWD.columns.values,
                           table4=[JPT.to_html(classes='table text-center bottommarg', index=False, escape=False)], title4=JPT.columns.values,
                           table5=[PUR.to_html(classes='table text-center bottommarg', index=False, escape=False)], title5=PUR.columns.values,
                           table6=[TSI.to_html(classes='table text-center bottommarg', index=False, escape=False)], title6=TSI.columns.values)
    else:
        # flash message if no results in dataframe rather than showing an empty one
        flash('No results found.')
        return redirect('/position')


@app.route('/Gene_name_results/<gene>', methods=['GET','POST'])
def gene_results(gene):
    data.gname=gene
    # query database for an gene name using pandas.query and store results in a new dataframe
    df2 = dff.query("Gene_Name == '%s'" %gene)
    # query database for a gene alias using pandas.query and store results in a new dataframe
    df2_2 = dff[dff['Aliases'].str.contains(gene)]
    # create list with both df's and concatanate them
    frames = [df2, df2_2]
    findf = pd.concat(frames)
    if findf.empty == False:
        # retrieve the positions of the first and last SNP's for a gene/alias and store in
        # the global Datastore class
        data.genepos1 = int(findf['POS'].iloc[0])
        data.genepos2 = int(findf['POS'].iloc[-1])
        # determine how many results we obtain
        hits = findf.shape[0]
        # hyperlink the ID of each row to their respective rsIDresults page
        findf['ID'] = findf['ID'].apply(lambda x: f'<a href="/rsIDresults/{x}">{x}</a>')
        # manipulate dataframe into multiple smaller df's for each population
        X = findf.iloc[:, 0:7]
        GIH = findf.iloc[:, [1, 2, 7, 8, 9, 10, 11]]
        GWD = findf.iloc[:, [1, 2, 12, 13, 14, 15, 16]]
        JPT = findf.iloc[:, [1, 2, 17, 18, 19, 20, 21]]
        PUR = findf.iloc[:, [1, 2, 22, 23, 24, 25, 26]]
        TSI = findf.iloc[:, [1, 2, 27, 28, 29, 30, 31]]
        # turning dataframes of each population into html tables using df.to_html
        return render_template('gene_results.html', hits=hits, table1=[X.to_html(classes='table text-center bottommarg', index=False, escape=False)], title1=X.columns.values,
                           table2=[GIH.to_html(classes='table text-center bottommarg', index=False, escape=False)], title2=GIH.columns.values,
                           table3=[GWD.to_html(classes='table text-center bottommarg', index=False, escape=False)], title3=GWD.columns.values,
                           table4=[JPT.to_html(classes='table text-center bottommarg', index=False,  escape=False)], title4=JPT.columns.values,
                           table5=[PUR.to_html(classes='table text-center bottommarg',  index=False, escape=False)], title5=PUR.columns.values,
                           table6=[TSI.to_html(classes='table text-center bottommarg', index=False, escape=False)], title6=TSI.columns.values)
    else:
        # flash message if no results in dataframe rather than showing an empty one
        flash('No results found.')
        return redirect('/Gene_name')



# read the genotype arrays as pickle files rather than dataframes
GIH = pd.read_pickle('GIH.pkl', compression={'method': 'zip'})
GWD = pd.read_pickle('GWD.pkl', compression={'method': 'zip'})
JPT = pd.read_pickle('JPT.pkl', compression={'method': 'zip'})
PUR = pd.read_pickle('PUR.pkl', compression={'method': 'zip'})
TSI = pd.read_pickle('TSI.pkl', compression={'method': 'zip'})

# defining the action of the stats route
@app.route('/stats', methods=['GET','POST'])
def stats():
    # assign the global variable set in /pos_results route into the stats function
    pos1 = data.pos1
    pos2 = data.pos2
    # get lists of chosen checkboxes
    req = request.form.getlist("pop")
    statreq = request.form.getlist("stat")
    # if/elif/else statement to send messages if errors are encountered
    if req ==[]:
        flash('Invalid selection(s), select at least one population and one summary statistic.', 'error')
        url = url_for('pos_results', position1=pos1, position2=pos2)
        url= url + '#Stats'
        return redirect(url)
    elif statreq == []:
        flash('Invalid selection(s), select at least one population and one summary statistic.', 'error')
        url = url_for('pos_results', position1=pos1, position2=pos2)
        url = url + '#Stats'
        return redirect(url)
    elif len(req) < 2 and 'd' in statreq:
        flash('Invalid selection, FST calculation requires at least 2 population selections.', 'error')
        url = url_for('pos_results', position1=pos1, position2=pos2)
        url = url + '#Stats'
        return redirect(url)
    else:
        # append pickle files to list if their respective populations are chosen
        pop_ls = []
        for i in req:
            if i == 'GIH':
                pop_ls.append(GIH)
            elif i == 'GWD':
                pop_ls.append(GWD)
            elif i == 'JPT':
                pop_ls.append(JPT)
            elif i == 'PUR':
                pop_ls.append(PUR)
            elif i== 'TSI':
                pop_ls.append(TSI)
        if request.method == 'POST':
            frames3 = []
            tabl2 = ''
            nucplot = ''
            homoplot=''
            tajplot=''
            fstplot = ''
            dd = []
            # receive input from form (for window size)
            wsize = request.form.get('window')
            #if statement to return message if condition not met
            if int(wsize) > (pos2-pos1):
                flash('Invalid window selection - outside range.', 'error')
                url = url_for('pos_results', position1=pos1, position2=pos2)
                url = url + '#Stats'
                return redirect(url)
            if int(wsize) == 0:
                flash('Invalid window selection - outside range.', 'error')
                url = url_for('pos_results', position1=pos1, position2=pos2)
                url = url + '#Stats'
                return redirect(url)
            #statistical functions are run if checkboxes for them are chosen
            if 'a' in statreq:
                aa = nucleotide_diversity(pop_ls, pos1, pos2)
                frames3.append(aa)
                nucplot = nucleotide_diversity_plot(pop_ls, pos1, pos2, int(wsize))
                #data2.nuc = nucleotide_diversity_plot(pop_ls, a, b, 1000)
            if 'b' in statreq:
                bb = homozygosity(pop_ls, pos1, pos2)
                frames3.append(bb)
                homoplot = homozygosity_plot(pop_ls, pos1, pos2, int(wsize))
                #data2.homo = homozygosity_plot(pop_ls, a, b, 1000)
            if 'c' in statreq:
                cc = tajima_D(pop_ls, pos1, pos2)
                frames3.append(cc)
                tajplot = tajima_D_plot(pop_ls, pos1, pos2, int(wsize))
                #data2.taj = tajima_D_plot(pop_ls, a, b, 1000)
            if 'd' in statreq:
                dd = fst(pop_ls, pos1, pos2)
                dd.index.rename('Population Code', inplace=True)
                tabl2 = dd.to_html(classes='table text-center bottommarg')
                fstplot = fst_plot(pop_ls, pos1, pos2, int(wsize))
                #data2.fst = fst_plot(pop_ls, a, b, 1000)
            # return FST results
            if frames3 == []:
                dd.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',')   # saves the fst results in a txt file
                return render_template('stats.html', tabl2=tabl2, plot4=fstplot)
            # statement to return  summary stats if fst is not chosen
            if len(dd) == 0:
                # concatanate dataframes in frames 3 list
                retdf = pd.concat(frames3, axis=1)
                retdf = retdf.fillna('-')
                retdf.index.rename('Population Code', inplace=True)
                # save dataframe as txt file
                retdf.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'),sep=',') # Save the data frame as a txt file into th Download folder
                return render_template('stats.html',tabl=[retdf.to_html(classes='table text-center bottommarg')],tit1=retdf.columns.values , plot1= nucplot, plot2=homoplot, plot3=tajplot)
            # return fst and summary stats (when fst is chosen)
            else:
                retdf = pd.concat(frames3, axis=1)
                retdf = retdf.fillna('-')
                retdf.index.rename('Population Code', inplace=True)
                retdf.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',') # Save the data frame as a txt file into th Download folder
                dd.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',', mode='a') # Join FSI data frame (dd) to the main dataframe (retdf) and save it as a txt file
                return render_template('stats.html',tabl=[retdf.to_html(classes='table text-center bottommarg')], tabl2=tabl2,tit1=retdf.columns.values , plot1= nucplot, plot2=homoplot, plot3=tajplot, plot4=fstplot)

        return render_template('stats.html')

@app.route('/statistics_for_gene', methods=['GET','POST'])
def genestats():
    # assign the global variable set in /gene_results route into the genestats function
    a = data.genepos1
    b = data.genepos2
    genename = str(data.gname)
    # get lists of chosen checkboxes
    req = request.form.getlist("pop")
    statreq = request.form.getlist("stat")
    # if/elif/else statement to send messages if errors are encountered
    if req == []:
        flash('Invalid selection(s), select at least one population and one summary statistic.', 'error')
        url = url_for('gene_results', gene=genename)
        url = url + '#Stats'
        return redirect(url)
    elif statreq == []:
        flash('Invalid selection(s), select at least one population and one summary statistic.', 'error')
        url = url_for('gene_results', gene=genename)
        url = url + '#Stats'
        return redirect(url)
    elif len(req) < 2 and 'd' in statreq:
        flash('Invalid selection, FST calculation requires at least 2 population selections.', 'error')
        url = url_for('gene_results', gene=genename)
        url = url + '#Stats'
        return redirect(url)
    else:
        # append pickle files to list if their respective populations are chosen
        pop_ls = []
        for i in req:
            if i == 'GIH':
                pop_ls.append(GIH)
            elif i == 'GWD':
                pop_ls.append(GWD)
            elif i == 'JPT':
                pop_ls.append(JPT)
            elif i == 'PUR':
                pop_ls.append(PUR)
            elif i == 'TSI':
                pop_ls.append(TSI)
        if request.method == 'POST':
            frames3 = []
            tabl2 = ''
            nucplot = ''
            homoplot = ''
            tajplot = ''
            fstplot = ''
            dd = []
            # receive input from form (for window size)
            wsize = request.form.get('window')
            # if statement to return message if condition not met
            if int(wsize) > (b - a):
                flash('Invalid window input - outside range.', 'error')
                url = url_for('gene_results', position1=a, position2=b)
                url = url + '#Stats'
                return redirect(url)
            if int(wsize) == 0:
                flash('Invalid window selection - outside range.', 'error')
                url = url_for('pos_results', position1=a, position2=b)
                url = url + '#Stats'
                return redirect(url)
            # statistical functions are run if checkboxes for them are chosen
            if 'a' in statreq:
                aa = nucleotide_diversity(pop_ls, a, b)
                frames3.append(aa)
                nucplot = nucleotide_diversity_plot(pop_ls, a, b, int(wsize))
            if 'b' in statreq:
                bb = homozygosity(pop_ls, a, b)
                frames3.append(bb)
                homoplot = homozygosity_plot(pop_ls, a, b, int(wsize))
            if 'c' in statreq:
                cc = tajima_D(pop_ls, a, b)
                frames3.append(cc)
                tajplot = tajima_D_plot(pop_ls, a, b, int(wsize))
            if 'd' in statreq:
                dd = fst(pop_ls, a, b)
                dd.index.rename('Population Code', inplace=True)
                tabl2 = dd.to_html(classes='table text-center bottommarg')
                fstplot = fst_plot(pop_ls, a, b, int(wsize))
            # return FST results
            if frames3 == []:
                dd.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',') # Save the data frame as a txt file into th Download folder
                return render_template('stats.html', tabl2=tabl2, plot4=fstplot)
            # statement to return  summary stats if fst is not chosen
            if len(dd) == 0:
                retdf = pd.concat(frames3, axis=1)
                retdf = retdf.fillna('-')
                retdf.index.rename('Population Code', inplace=True)
                retdf.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'),sep=',') # Save the data frame as a txt file into th Download folder
                # save dataframe as txt
                return render_template('stats.html',tabl=[retdf.to_html(classes='table text-center bottommarg')],tit1=retdf.columns.values , plot1= nucplot, plot2=homoplot, plot3=tajplot)
            # return fst and summary stats (when fst is chosen)
            else:
                retdf = pd.concat(frames3, axis=1)
                retdf = retdf.fillna('-')
                retdf.index.rename('Population Code', inplace=True)
                retdf.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',') # Save the data frame as a txt file into th Download folder
                dd.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',', mode='a') # Join FSI data frame (dd) to the main dataframe (retdf) and save it as a txt file
                return render_template('stats.html',tabl=[retdf.to_html(classes='table text-center bottommarg')], tabl2=tabl2,tit1=retdf.columns.values , plot1= nucplot, plot2=homoplot, plot3=tajplot, plot4=fstplot)

        return render_template('stats.html')


@app.route('/return-files/')
def return_files():
    """ Download the Results.txt file and delete it afterwards """
    try:
        # Download the Results.txt file
        return send_file(app.config['path_Results.txt'],  as_attachment= True), \
               os.remove(os.path.join(app.config['path_Results.txt']))  # Once the file is downloaded delete it

    except Exception:
        flash('Download error, please try again.', 'error')
        return  "Unable to download file, please try again." # if problems arise display this string in new window

@scheduler.task('interval', id='del_files', seconds=1800, misfire_grace_time=None)
def del_files():
    """ Delete the saved statistics results txt files after 30 minutes """
    current_time = time.time() # Take the current time
    size_Download_folder = os.stat(app.config['path_Download_folder']).st_size # Check the size of the download folder

    if size_Download_folder > 0: # if the download folder is not empty

        for stat_results_filename in os.listdir(app.config['path_Download_folder']):  # For each file in the download folder, although there will only be one file.
            if os.path.getmtime(os.path.join(app.config['path_Download_folder'], stat_results_filename)) < current_time  -1:  # Check if they are less than  1 day
                if os.path.isfile(os.path.join(app.config['path_Download_folder'], stat_results_filename)):  # check if the file exists in the file system
                    os.remove(os.path.join(app.config['path_Download_folder'], stat_results_filename)) # Remove the file
    else:
        pass


# start the web server
if __name__ == '__main__':
    scheduler.start() # Start app schedular once the app has been started
    atexit.register(lambda: scheduler.shutdown()) # Shut down apscheduler after leaving the application
    app.run(debug=True)
