 ############################## imports #################################
import os
from flask_apscheduler import APScheduler
import atexit
import time

 ############################## Set up #################################
 
# Set the path to the download folder
app.config['path_Download_folder'] = os.path.join(os.getcwd(), 'Downloads')

# Set configuration for app schedualer
class Config:
    SCHEDULER_API_ENABLED = True

# Initialize scheduler
scheduler = APScheduler()


############################ Function  ###################################

@scheduler.task('interval', id='del_files', seconds=1800, misfire_grace_time=None)
def del_files():
    """ Delete the saved statistics  results  csv files after 30 minutes """
    current_time = time.time()
    # Check the size of the download folder 
    size_Download_folder = os.stat(app.config['path_Download_folder']).st_size

    if size_Download_folder > 0:
        # For each file in the download folder, although there will only be one file. This was done thinking in a future inplementation.
        for stat_results_filename in os.listdir(app.config['path_Download_folder']):
            # Check if they are less than  1 day
            if os.path.getmtime(os.path.join(app.config['path_Download_folder'], stat_results_filename)) < current_time  -1:
                # check if the file exists in the file system
                if os.path.isfile(os.path.join(app.config['path_Download_folder'], stat_results_filename)):
                    # Remove the file 
                    os.remove(os.path.join(app.config['path_Download_folder'], stat_results_filename))
    else:
        pass

      ############################ Extra code in the start app code ###################################
      
      if __name__ == '__main__':
        # Start app schedualer once the app has been started
        scheduler.start()
        # Shut down apscheduler after leaving the application
        atexit.register(lambda: scheduler.shutdown())
        app.run(debug=True)
      
      
      
      
