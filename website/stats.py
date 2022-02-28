import allel
import numpy as np
import pandas as pd
from allel.stats.window import moving_statistic
import plotly.graph_objects as go
import plotly as py
import plotly
import pandas as pd
import plotly.express as px
import plotly.io as pio
import json
import os

###### Nucleotide Diversity (accepts a list of pkl) ######

def nucleotide_diversity(filelist, mini, maxi):
    """ Computes nucleotide diversity within a given region from a list of objects of pkl files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = file
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]  # collapse dataframe for region required
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]  # extract the population code
        arr = pd.DataFrame(df).to_numpy()  # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)
        geno = haplotypes.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        pi = allel.sequence_diversity(position, ac, start=mini, stop=maxi)  # nucleotide diversity calculation
        results[country] = pi
    results = pd.DataFrame.from_dict(results, orient='index',
                                     columns=['Nucleotide Diversity'])  # return all results as a dataframe
    return results


###### Tajima's D (accepts a list of pkl)######

def tajima_D(filelist, mini, maxi):
    """ Computes Tajima's D within a given region from list of objects of pkl files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = file
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]  # collapse dataframe for region required
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]  # extract the population code
        arr = pd.DataFrame(df).to_numpy()  # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)
        geno = haplotypes.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        t = allel.tajima_d(ac, pos=position, start=mini, stop=maxi)  # Tajima's D calculation
        results[country] = t
    results = pd.DataFrame.from_dict(results, orient='index', columns=['Tajima D'])  # return all results as a dataframe
    return results

    ###### Homozygosity (accepts a list of pkl) ######


def homozygosity(filelist, mini, maxi):
    """ Takes a list of objects of pkl files as its input and returns the average homozygosity over a given region. """
    results = {}
    for file in filelist:
        df = file
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]  # collapse data for region required
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]
        num = df.shape[1] / 2
        arr = pd.DataFrame(df).to_numpy()
        haplotypes = allel.HaplotypeArray(arr)  # haplotype array
        gt = haplotypes.to_genotypes(ploidy=2)  # convert to genotype array
        h = gt.is_hom()
        homo = np.sum(h, axis=1)  # homozygous counts
        homo = pd.DataFrame(homo, columns=['HOMOZYGOSITY'])
        homo['HOMOZYGOSITY'] = homo['HOMOZYGOSITY'].div(num)  # homozygous frequency calculation
        homozygosity = homo['HOMOZYGOSITY'].sum() / len(homo.index)  # average homozygosity across region
        results[country] = homozygosity
    results = pd.DataFrame.from_dict(results, orient='index', columns=['Homozygosity'])  # return results as a dataframe
    return results


###### FST (accepts a list of pkl) ######

def fst(filelist, mini, maxi):
    pairs = [(filelist[i], filelist[j]) for i in range(len(filelist)) for j in
             range(i + 1, len(filelist))]  ##get all the possible combinations if the user pick more than 2 populations
    inputlist1 = [item[0] for item in pairs]  ### seperate the combination lists into 2 list
    inputlist2 = [item[1] for item in pairs]
    col = []
    col1 = []
    # newData = [tuple(map(lambda i: str.replace(i, ".pkl"," "), tup)) for tup in pairs]
    output = []  ## create an empty list and append the results to the list to avoid the return statament to break the for loop
    for (x, y) in zip(inputlist1, inputlist2):  ## a for loop to compute Allele frequency of each element of the lists
        df = x
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        arr = pd.DataFrame(df).to_numpy()  # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)
        geno = haplotypes.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        results = df.columns[0]
        col.append(results)
        df1 = y
        df1 = df1.loc[(df1['POS'] >= mini) & (df1['POS'] <= maxi)]
        new_dataframe1 = df1.filter(['POS'])  # extract the position from the dataframe
        position1 = pd.DataFrame(new_dataframe1, columns=['POS']).to_numpy()
        position1 = position1.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df1 = df1.drop(['POS'], axis=1)
        arr1 = pd.DataFrame(df1).to_numpy()  # haplotype dataframe to numpy array
        haplotypes1 = allel.HaplotypeArray(arr1)
        geno1 = haplotypes1.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac1 = geno1.count_alleles()
        results1 = df1.columns[0]
        col1.append(results1)
        num, dem = allel.hudson_fst(ac, ac1)  ## compute fst using allel package
        fst = np.sum(num) / np.sum(dem)
        fst_lst = fst.tolist()  ## tranform output array to a list in order to append it to the output list
        output.append(fst_lst)
        df_fst = pd.DataFrame(output)  ## trasnfom the output list to a data frame for easier reading
        combi = list(zip(col, col1))

    fst_ = df_fst.T
    combi = [map(str, c) for c in combi]
    combi = [(' vs '.join(c)) for c in combi]
    fst_.columns = combi
    fst_ = fst_.rename(columns=lambda x: x.strip('*'))
    fst_final = fst_.T
    fst_final.columns = ['Fixation Index']
    return fst_final

    ###### Nucleotide Diversity PLOT ######


def nucleotide_diversity_plot(filelist, mini, maxi, size):
    """ Computes nucleotide diversity within a given region and a given window size  from a list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = (file)
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]  # extract the population code
        arr = pd.DataFrame(df).to_numpy()  # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)
        geno = haplotypes.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        pi, windows, n_bases, counts = allel.windowed_diversity(position, ac, size, start=mini,
                                                                stop=maxi)  # nucleotide diversity calculation
        results[country] = pi
    results = pd.DataFrame.from_dict(results, orient='index').T  # return all results as a dataframe
    windows = windows.tolist()  ##transform windows array to a list
    windowsdf = pd.DataFrame(windows)
    final_df = pd.concat([windowsdf, results],
                         axis=1)  ##create a data frame with the windows and the nucleotide diversity scores
    windows_axis = (final_df.iloc[:, 2:])  ## set the windows position as the the x-axis
    ## plot nucleotide diversity scores by windows
    fig = go.Figure()
    for idx, col in enumerate(windows_axis.columns, 0):
        fig.add_trace(go.Scatter(x=final_df.iloc[:, 0], y=(windows_axis.iloc[:, idx]), mode='lines', name=col))
        fig.update_layout(
            title_text="Nucleotide Diversity across windows")
        fig.update_layout(
            xaxis=dict(title="Window Positions (bp)"),
            yaxis=dict(title="Nucleotide Diversity"))

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON


###### Homozygosity PLOT ######

def homozygosity_plot(filelist, mini, maxi, size):
    """ Takes a CSV file as its input and returns the average homozygosity over a given region. """
    results = {}
    for file in filelist:
        df = (file)
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]
        arr = pd.DataFrame(df).to_numpy()  ## convert the dataframe into a numpy array
        haplotypes = allel.HaplotypeArray(arr)  ## extract the haplotype array
        genotype = haplotypes.to_genotypes(ploidy=2)
        num = genotype.n_samples  ## extract the number of samples by populations
        homozygote_count = np.sum(genotype.is_hom(), axis=1)  # total number of homozygote
        homozygote_count = pd.DataFrame(homozygote_count,
                                        columns=['HOMOZYGOSITY'])  ## convert the array into a dataframe
        homozygote_count['HOMOZYGOSITY'] = homozygote_count['HOMOZYGOSITY'].div(
            num)  ## divide homozygote by the num to get frequency of the homozygote genotype
        homozygote_count_array = homozygote_count.to_numpy()
        a, window, count = allel.windowed_statistic(position, homozygote_count_array, np.mean, size, start=mini,
                                                    stop=maxi)  ##compute the mean of homozygote indivduals in windows
        results[country] = a
    results = pd.DataFrame.from_dict(results,
                                     orient='index').T  # return all results as a dataframe and transpose it for easy manipulations
    windows_lst = window.tolist()  ##transform windows array to a list
    windowsdf = pd.DataFrame(windows_lst)
    final_df = pd.concat([windowsdf, results],
                         axis=1)  ##create a data frame with the windows and the nucleotide diversity scores
    final_df.fillna(0)
    windows_axis = (final_df.iloc[:, 2:])  ## set the windows position as the the x-axis
    ## plot nucleotide diversity scores by windows
    fig = go.Figure()
    for idx, col in enumerate(windows_axis.columns, 0):
        fig.add_trace(go.Scatter(x=final_df.iloc[:, 0], y=windows_axis.iloc[:, idx], mode='lines', name=col))
        fig.update_layout(
            title_text="Homozygosity across windows")
        fig.update_layout(
            xaxis=dict(title="Window Positions (bp)"),
            yaxis=dict(title="Homozygosity"))

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON


###### Tajima D PLOT ######

def tajima_D_plot(filelist, mini, maxi, size):
    """ Computes Tajima's D within a given region from list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = file
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]  # extract the population code
        arr = pd.DataFrame(df).to_numpy()  # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)
        geno = haplotypes.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        D, windows, counts = allel.windowed_tajima_d(position, ac, size=size, start=mini, stop=maxi,
                                                     min_sites=0)  # Tajima's D calculation
        results[country] = D
    results = pd.DataFrame.from_dict(results,
                                     orient='index').T  # return all results as a dataframe and transpose it for easy manipulations
    results_final = results.fillna(0)  ## fill the Nan value with 0
    # windows=windows.tolist() ##transform windows array to a list
    windowsdf = pd.DataFrame(windows)
    final_df = pd.concat([windowsdf, results_final],
                         axis=1)  ##create a data frame with the windows and the nucleotide diversity scores
    windows_axis = (final_df.iloc[:, 2:])  ## set the windows position as the the x-axis
    ## plot nucleotide diversity scores by windows
    fig = go.Figure()
    for idx, col in enumerate(windows_axis.columns, 0):
        fig.add_trace(go.Scatter(x=final_df.iloc[:, 0], y=windows_axis.iloc[:, idx], mode='lines', name=col))
        fig.update_layout(
            title_text="TajimaD  across windows")
        fig.update_layout(
            xaxis=dict(title="Window Positions (bp)"),
            yaxis=dict(title="TajimaD"))

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON


###### FST PLOT ######

def fst_plot(filelist, mini, maxi, size):
    pairs = [(filelist[i], filelist[j]) for i in range(len(filelist)) for j in
             range(i + 1, len(filelist))]  ##get all the possible combinations if the user pick more than 2 populations
    inputlist1 = [item[0] for item in pairs]  ### seperate the combination lists into 2 list
    inputlist2 = [item[1] for item in pairs]
    output = []
    col = []
    col1 = []
    for (x, y) in zip(inputlist1, inputlist2):  ## a for loop to compute Allele frequency of each element of the lists
        df = (x)
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        arr = pd.DataFrame(df).to_numpy()  # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)
        geno = haplotypes.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        results = df.columns[0]
        col.append(results)
        df1 = (y)
        df1 = df1.loc[(df1['POS'] >= mini) & (df1['POS'] <= maxi)]
        new_dataframe1 = df1.filter(['POS'])  # extract the position from the dataframe
        position1 = pd.DataFrame(new_dataframe1, columns=['POS']).to_numpy()
        position1 = position1.flatten()  # alter position to one dimensional numpy array for correct input for sequence diversity function
        df1 = df1.drop(['POS'], axis=1)
        arr1 = pd.DataFrame(df1).to_numpy()  # haplotype dataframe to numpy array
        haplotypes1 = allel.HaplotypeArray(arr1)
        geno1 = haplotypes1.to_genotypes(
            ploidy=2)  # reshape haplotype array to genotype array for subsequent correct function input
        ac1 = geno1.count_alleles()
        results1 = df1.columns[0]
        col1.append(results1)
        num, den = allel.hudson_fst(ac, ac1)  ## compute fst using allel package

        def average_fst(wn, wd):
            return np.nansum(wn) / np.nansum(wd)

        # calculate average Fst in windows
        fst, windows, counts = allel.windowed_statistic(position, values=(num, den), statistic=average_fst, size=size,
                                                        start=mini, stop=maxi)
        fst_lst = fst.tolist()  ## tranform output array to a list in order to append it to the output list
        output.append(fst_lst)
        windows_lst = windows.tolist()
        # fst_lst=fst.tolist() ## tranform output array to a list in order to append it to the output list
        # output.append(fst_lst)
        df_fst = pd.DataFrame(output)  ## trasnfom the output list to a data frame for easier reading
        combi = list(zip(col, col1))
    fst_ = df_fst.T
    combi = [map(str, c) for c in combi]
    combi = [(' vs '.join(c)) for c in combi]
    fst_.fillna(0)
    fst_.columns = combi

    windowsdf = pd.DataFrame(windows_lst)
    a = windowsdf.iloc[:, 0]

    fig = go.Figure()
    for idx, col in enumerate(fst_.columns, 0):
        fig.add_trace(go.Scatter(x=a, y=(fst_.iloc[:, idx]), mode='lines', name=str(col)))
        fig.update_layout(title_text="FST")
        fig.update_layout(xaxis=dict(title="Window Positions (bp)"),
                          yaxis=dict(title="Fixation Index(FST)"))

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON