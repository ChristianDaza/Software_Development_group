import allel
import numpy as np
import modin.pandas as pd

###### Nucleotide Diversity ######

def nucleotide_diversity(filelist, mini, maxi):
    """ Computes nucleotide diversity within a given region from a list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = pd.read_csv(file, index_col=0)
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]    # extract the population code
        arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)   
        geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        pi = allel.sequence_diversity(position, ac, start=mini, stop=maxi) # nucleotide diversity calculation
        results[country] = pi
    results = pf.DataFrame.from_dict(results, orient='index', columns=['Nucleotide Diversity'])  # return all results as a dataframe
    return results
  
  
###### Haplotype Diversity ######

def haplotype_diversity(filelist, mini, maxi):
    """ Computes Haplotype Diversity within a given region from a list of CSV files containing position and haplotype data. """
    results ={}
    for file in filelist:
        df = pd.read_csv(file, index_col=0)
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)] # retain haplotype data for only within the region specified
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]    # extract the population code
        arr = pd.DataFrame(df).to_numpy()   # reshape dataframe to a numpy array
        haplotypes = allel.HaplotypeArray(arr) # numpy array to haplotype array
        hd = allel.haplotype_diversity(haplotypes)  # haplotype diversity calculation
        results[country]= hd
    results = pd.Dataframe.from_dict(results, orient='index', columns=['Haplotype Diversity'])  # return all results as a dataframe
    return results 
  

###### Tajima's D ######

def tajima_D(filelist, mini, maxi):
    """ Computes Tajima's D within a given region from list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = pd.read_csv(file, index_col=0)
        df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]    # extract the population code
        arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
        haplotypes = allel.HaplotypeArray(arr)   
        geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
        ac = geno.count_alleles()
        t = allel.tajima_d(ac, pos=position, start=mini, stop=maxi)  # Tajima's D calculation
        results[country] = t
    results = pd.Dataframe.from_dict(results, orient='index', columns=['Tajima D'])  # return all results as a dataframe
    return results

###### Homozygosity ######

def homozygosity(inputfile, mini, maxi):
    """ Takes a CSV file as its input and returns the average homozygosity over a given region. """
    df = pd.read_csv(inputfile, index_col=0)
    df= df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
    new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
    position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
    position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
    df = df.drop(['POS'], axis=1)
    num = df.shape[1] / 2
    arr = pd.DataFrame(df).to_numpy()
    haplotypes = allel.HaplotypeArray(arr)
    gt = haplotypes.to_genotypes(ploidy=2)
    h = gt.is_hom()
    homo = np.sum(h, axis=1)
    homo= pd.DataFrame(homo, columns=['HOMOZYGOSITY'])
    homo['HOMOZYGOSITY'] = homo['HOMOZYGOSITY'].div(num)
    homozygosity=  homo['HOMOZYGOSITY'].sum() / len(homo.index)
    return homozygosity
