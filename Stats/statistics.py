import allel
import numpy as np
import modin.pandas as pd

###### nucleotide diversity ######

def nucleotide_diversity(inputfile, mini, maxi):
    """ Computes nucleotide diversity within a given region from a CSV file containing position and haplotype data. """
    df = pd.read_csv(inputfile, index_col=0)
    new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
    position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
    position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
    df = df.drop(['POS'], axis=1)
    arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
    haplotypes = allel.HaplotypeArray(arr)   
    geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
    ac = geno.count_alleles()
    pi = allel.sequence_diversity(position, ac, start=mini, stop=maxi) # nucleotide diversity calculation
    return pi
  
  
###### Haplotype Diversity ######

def haplotype_diversity(inputfile, mini, maxi):
    """ Computes Haplotype Diversity within a given region from a CSV file containing position and haplotype data. """
    df = pd.read_csv(inputfile, index_col=0)
    df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)] # retain haplotype data for only within the region specified
    df = df.drop(['POS'], axis=1)
    arr = pd.DataFrame(df).to_numpy()   # reshape dataframe to a numpy array
    haplotypes = allel.HaplotypeArray(arr) # numpy array to haplotype array
    hd = allel.haplotype_diversity(haplotypes)  # haplotype diversity calculation
    return hd 
  

###### Tajima's D ######

def tajima_D(inputfile, mini, maxi):
    """ Computes Tajima's D within a given region from a CSV file containing position and haplotype data. """
    df = pd.read_csv(inputfile, index_col=0)
    df = df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
    new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
    position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
    position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
    df = df.drop(['POS'], axis=1)
    arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
    haplotypes = allel.HaplotypeArray(arr)   
    geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
    ac = geno.count_alleles()
    t = allel.tajima_d(ac, pos=position, start=mini, stop=maxi)  # Tajima's D calculation
    return t
