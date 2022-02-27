########################################### Imports ################################################
import pandas as pd
import os

########################################### Setting up paths #######################################

# Set the path for the  saved  statistical results file
app.config['path_Results.txt'] = os.path.join(os.getcwd(), 'Downloads', 'Results.txt')

###################### Saving the data frames  into  the Download folder ###########################

# Save the data frame as a csv file into th Download folder
retdf.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',')

# Join FSI data frame (dd) to the main dataframe (retdf)
 dd.to_csv(path_or_buf= os.path.join(os.getcwd(), 'Downloads', 'Results.txt'), sep=',', mode='a')
 

########################################### Download fucntion #######################################

@app.route("/return_files")
def return_files():
    """ Download the Results.txt file and delet it afterwards """
    try:
         # Download the Results.txt file
        return send_file(app.config['path_Results.txt'],  as_attachment= True), \
              # Once the file is downloaded delete it
	      os.remove(os.path.join(app.config['path_Results.txt'])),

    except FileNotFoundError:
        return "Unable to download file please try again"


