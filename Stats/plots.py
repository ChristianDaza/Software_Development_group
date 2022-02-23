import allel
import numpy as np
import pandas as pd
from allel.stats.window import moving_statistic
import plotly.graph_objects as go
import plotly as py
import pandas as pd
import plotly.express as px
import plotly.io as pio
size=1000 ## a default window size 
def nucleotide_diversity_by_window(filelist,size, mini, maxi):
    """ Computes nucleotide diversity within a given region and a given window size  from a list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = pd.read_pickle(file)
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]    # extract the population code
        arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)   
        geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        pi, windows, n_bases, counts = allel.windowed_diversity(position, ac,size, start=mini, stop=maxi)# nucleotide diversity calculation
        results[country] = pi
    results = pd.DataFrame.from_dict(results, orient='index').T  # return all results as a dataframe
    windows=windows.tolist() ##transform windows array to a list 
    windowsdf=pd.DataFrame(windows)    
    final_df=pd.concat([windowsdf,results], axis=1)##create a data frame with the windows and the nucleotide diversity scores 
    windows_axis=(final_df.iloc[: , 2:]) ## set the windows position as the the x-axis
    ## plot nucleotide diversity scores by windows 
    fig = go.Figure()
    for idx, col in enumerate(windows_axis.columns, 0):
        fig.add_trace(go.Scatter(x = final_df.iloc[: , 0] , y = np.log(windows_axis.iloc[:,idx]), mode ='lines', name = 'nd'+col))
        fig.update_layout(
    title_text="Nucletotide Diversity across windows")
        fig.update_layout(
    xaxis=dict(title="Window Positions (bp)"),
        yaxis=dict(title="Nucletotide Diversity"))
    Nucleotide_diversity_fig=fig.show()

    return  Nucleotide_diversity_fig
