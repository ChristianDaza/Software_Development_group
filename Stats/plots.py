import allel
import numpy as np
import pandas as pd
from allel.stats.window import moving_statistic
import plotly.graph_objects as go
import plotly as py
import pandas as pd
import plotly.express as px
import plotly.io as pio




def fst(filelist,mini,maxi):
        pairs=[(filelist[i],filelist[j]) for i in range(len(filelist)) for j in range(i+1, len(filelist))]##get all the possible combinations if the user pick more than 2 populations 
        inputlist1=[item[0] for item in pairs] ### seperate the combination lists into 2 list 
        inputlist2=[item[1] for item in pairs]
        col=[]
        col1=[]
        #newData = [tuple(map(lambda i: str.replace(i, ".pkl"," "), tup)) for tup in pairs]
        output=[] ## create an empty list and append the results to the list to avoid the return statament to break the for loop 
        for (x,y) in zip(inputlist1,inputlist2):## a for loop to compute Allele frequency of each element of the lists 
                df = x#, index_col=0)
                df = df.loc[(df['POS'] >=  mini ) & (df['POS'] <= maxi)]
                new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
                position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
                position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
                df = df.drop(['POS'], axis=1)
                arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
                haplotypes = allel.HaplotypeArray(arr)   
                geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
                ac= geno.count_alleles()
                results=df.columns[0]
                col.append(results)
                df1 = y#, index_col=0)
                df1 = df1.loc[(df1['POS'] >=  mini ) & (df1['POS'] <= maxi)]
                new_dataframe1 = df1.filter(['POS'])  # extract the position from the dataframe
                position1 = pd.DataFrame(new_dataframe1, columns=['POS']).to_numpy()
                position1 = position1.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
                df1 = df1.drop(['POS'], axis=1)
                arr1 = pd.DataFrame(df1).to_numpy()     # haplotype dataframe to numpy array
                haplotypes1 = allel.HaplotypeArray(arr1)   
                geno1 = haplotypes1.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
                ac1= geno1.count_alleles()
                results1=df1.columns[0]
                col1.append(results1)
                num,dem= allel.hudson_fst(ac,ac1) ## compute fst using allel package 
                fst = np.sum(num) / np.sum(dem)
                fst_lst=fst.tolist() ## tranform output array to a list in order to append it to the output list 
                output.append(fst_lst)
                df_fst = pd.DataFrame (output) ## trasnfom the output list to a data frame for easier reading 
                combi=list(zip(col ,col1))
        
        fst_=df_fst.T
        fst_.columns=combi
        fst_final=fst_.T
        fst_final.columns = ['Fixation Index']
        return fst_final
        



            


size=1000 
def fst_window_plot(filelist,size,mini,maxi):
    
        pairs=[(filelist[i],filelist[j]) for i in range(len(filelist)) for j in range(i+1, len(filelist))]##get all the possible combinations if the user pick more than 2 populations 
        inputlist1=[item[0] for item in pairs] ### seperate the combination lists into 2 list 
        inputlist2=[item[1] for item in pairs]
        output=[] 
        col=[]
        col1=[]
        for (x,y) in zip(inputlist1,inputlist2): ## a for loop to compute Allele frequency of each element of the lists 
                df = (x)
                df = df.loc[(df['POS'] >=  mini ) & (df['POS'] <= maxi)]
                new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
                position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
                position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
                df = df.drop(['POS'], axis=1)
                arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
                haplotypes = allel.HaplotypeArray(arr)   
                geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
                ac= geno.count_alleles()
                results=df.columns[0]
                col.append(results)
                df1 = (y)
                df1 = df1.loc[(df1['POS'] >=  mini ) & (df1['POS'] <= maxi)]
                new_dataframe1 = df1.filter(['POS'])  # extract the position from the dataframe
                position1 = pd.DataFrame(new_dataframe1, columns=['POS']).to_numpy()
                position1 = position1.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
                df1 = df1.drop(['POS'], axis=1)
                arr1 = pd.DataFrame(df1).to_numpy()     # haplotype dataframe to numpy array
                haplotypes1 = allel.HaplotypeArray(arr1)   
                geno1 = haplotypes1.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
                ac1= geno1.count_alleles()
                results1=df1.columns[0]
                col1.append(results1)
                num,den= allel.hudson_fst(ac,ac1) ## compute fst using allel package 
                def average_fst(wn, wd):
                    return np.nansum(wn) / np.nansum(wd)
                # calculate average Fst in windows
                fst, windows, counts = windowed_statistic(position, values=(num, den),statistic=average_fst,size=size, start=mini,stop=maxi)
                fst_lst=fst.tolist() ## tranform output array to a list in order to append it to the output list 
                output.append(fst_lst)
                windows_lst=windows.tolist()
                #fst_lst=fst.tolist() ## tranform output array to a list in order to append it to the output list 
                #output.append(fst_lst)
                df_fst = pd.DataFrame (output)## trasnfom the output list to a data frame for easier reading 
                combi=list(zip(col ,col1))
        fst_=df_fst.T
        fst_.columns=combi
        #fst_final=fst_.T
        #fst_final.columns = ['Fixation Index']
    
                
        windowsdf=pd.DataFrame(windows_lst)
        a=windowsdf.iloc[:, 0]

        #fin.columns = 
        fig = go.Figure()
        for idx, col in enumerate(fst_.columns, 0):
            fig.add_trace(go.Scatter(x = a , y = (fst_.iloc[:,idx]), mode ='lines', name = str(col)))
            fig.update_layout(title_text="FST")
            fig.update_layout(xaxis=dict(title="Window Positions (bp)"),
        yaxis=dict(title="Fixation Index(FST)"))
    
        return fig.show()
  
         
    
            



  
         
    
            



size=1000 ## a default window size 

size=1000 ## a default window size 
def nucleotide_diversity_by_window(filelist,size, mini, maxi):
    """ Computes nucleotide diversity within a given region and a given window size  from a list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = (file)
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
        fig.add_trace(go.Scatter(x = final_df.iloc[: , 0] , y = np.log(windows_axis.iloc[:,idx]), mode ='lines', name = col))
        fig.update_layout(
    title_text="Nucletotide Diversity across windows")
        fig.update_layout(
    xaxis=dict(title="Window Positions (bp)"),
        yaxis=dict(title="Nucletotide Diversity"))
    Nucleotide_diversity_fig=fig.show()

    return  Nucleotide_diversity_fig








def homozygosit_by_window(filelist,size, mini, maxi):
    """ Takes a CSV file as its input and returns the average homozygosity over a given region. """
    results = {}
    for file in filelist:
        df = pd.read_pickle(file)
        df= df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]
        arr = pd.DataFrame(df).to_numpy() ## convert the dataframe into a numpy array 
        haplotypes = allel.HaplotypeArray(arr) ## extract the haplotype array 
        genotype = haplotypes.to_genotypes(ploidy=2)
        num=genotype.n_samples## extract the number of samples by populations 
        homozygote_count = np.sum(genotype.is_hom(), axis=1) #total number of homozygote 
        homozygote_count= pd.DataFrame(homozygote_count, columns=['HOMOZYGOSITY']) ## convert the array into a dataframe 
        homozygote_count['HOMOZYGOSITY'] = homozygote_count['HOMOZYGOSITY'].div(num)## divide homozygote by the num to get frequency of the homozygote genotype 
        homozygote_count_array=homozygote_count.to_numpy()
        a,window,count=allel.windowed_statistic(position,homozygote_count_array,np.mean,size, start=mini, stop=maxi) ##compute the mean of homozygote indivduals in windows 
        results[country] = a
        results = pd.DataFrame.from_dict(results, orient='index').T # return all results as a dataframe and transpose it for easy manipulations 
        windows_lst=window.tolist() ##transform windows array to a list 
        windowsdf=pd.DataFrame(windows_lst)    
        final_df=pd.concat([windowsdf,results], axis=1)##create a data frame with the windows and the nucleotide diversity scores 
        windows_axis=(final_df.iloc[: , 2:]) ## set the windows position as the the x-axis
    ## plot nucleotide diversity scores by windows 
    fig = go.Figure()
    for idx, col in enumerate(windows_axis.columns, 0):
        fig.add_trace(go.Scatter(x = final_df.iloc[: , 0] , y = windows_axis.iloc[:,idx], mode ='lines', name = col))
        fig.update_layout(
    title_text="Homozygosity across windows")
        fig.update_layout(
    xaxis=dict(title="Window Positions (bp)"),
        yaxis=dict(title="Homozygosity"))
    Homozygosity_fig=fig.show()
    return Homozygosity_fig




size=1000 ## a default window size 
def tajima_D_by_window1(filelist,size, mini, maxi):
    """ Computes nucleotide diversity within a given region and a given window size  from a list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = (file)
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]    # extract the population code
        arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)   
        geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        pi, windows, counts = allel.windowed_tajima_d(position, ac,size, start=mini, stop=maxi,min_sites=0)# nucleotide diversity calculation
        results[country] = pi
    results = pd.DataFrame.from_dict(results, orient='index').T  # return all results as a dataframe
    windows=windows.tolist() ##transform windows array to a list 
    windowsdf=pd.DataFrame(windows)    
    final_df=pd.concat([windowsdf,results], axis=1)##create a data frame with the windows and the nucleotide diversity scores 
    windows_axis=(final_df.iloc[: , 2:]) ## set the windows position as the the x-axis

    fig = go.Figure()
    for idx, col in enumerate(windows_axis.columns, 0):
        fig.add_trace(go.Scatter(x = final_df.iloc[: , 0] , y = (windows_axis.iloc[:,idx]), mode ='lines', name = col))
        fig.update_layout(
    title_text="TajimaD across windows")
        fig.update_layout(
    xaxis=dict(title="Window Positions (bp)"),
        yaxis=dict(title="TajimaD"))
    Nucleotide_diversity_fig=fig.show()

    return  Nucleotide_diversity_fig
