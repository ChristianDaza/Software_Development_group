    
###### Nucleotide Diversity (accepts a list of pkl) ######

def nucleotide_diversity_pkl(filelist, mini, maxi):
    """ Computes nucleotide diversity within a given region from a list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = file
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
    results = pd.DataFrame.from_dict(results, orient='index', columns=['Nucleotide Diversity'])  # return all results as a dataframe
    return results
  
 ###### Tajima's D (accepts a list of pkl)######

def tajima_D_pkl(filelist, mini, maxi):
    """ Computes Tajima's D within a given region from list of CSV files containing position and haplotype data. """
    results = {}
    for file in filelist:
        df = file
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
    results = pd.DataFrame.from_dict(results, orient='index', columns=['Tajima D'])  # return all results as a dataframe
    return results
  
  ###### Homozygosity (accepts a list of pkl) ######

  def homozygosity_pkl(filelist, mini, maxi):
    """ Takes a list of CSV files as its input and returns the average homozygosity over a given region. """
    results = {}
    for file in filelist:
        df = file
        df= df.loc[(df['POS'] >= mini) & (df['POS'] <= maxi)]
        new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
        position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
        position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
        df = df.drop(['POS'], axis=1)
        country = df.columns[0]
        num = df.shape[1] / 2
        arr = pd.DataFrame(df).to_numpy()
        haplotypes = allel.HaplotypeArray(arr)
        gt = haplotypes.to_genotypes(ploidy=2)
        h = gt.is_hom()
        homo = np.sum(h, axis=1)
        homo= pd.DataFrame(homo, columns=['HOMOZYGOSITY'])
        homo['HOMOZYGOSITY'] = homo['HOMOZYGOSITY'].div(num)
        homozygosity =  homo['HOMOZYGOSITY'].sum() / len(homo.index)
        results[country] = homozygosity
    results = pd.DataFrame.from_dict(results, orient='index', columns=['Homozygosity'])
    return results
   
 ###### FST (accepts a list of pkl) ######

  def fst(filelist,mini,maxi):
        pairs=[(filelist[i],filelist[j]) for i in range(len(filelist)) for j in range(i+1, len(filelist))]##get all the possible combinations if the user pick more than 2 populations 
        inputlist1=[item[0] for item in pairs] ### seperate the combination lists into 2 list 
        inputlist2=[item[1] for item in pairs]
        newData = [tuple(map(lambda i: str.replace(i, ".pkl"," "), tup)) for tup in pairs]
        output=[] ## create an empty list and append the results to the list to avoid the return statament to break the for loop 
        for (x,y) in zip(inputlist1,inputlist2):## a for loop to compute Allele frequency of each element of the lists 
                df = pd.read_pickle(x)#, index_col=0)
                df = df.loc[(df['POS'] >=  mini ) & (df['POS'] <= maxi)]
                new_dataframe = df.filter(['POS'])  # extract the position from the dataframe
                position = pd.DataFrame(new_dataframe, columns=['POS']).to_numpy()
                position = position.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
                df = df.drop(['POS'], axis=1)
                arr = pd.DataFrame(df).to_numpy()     # haplotype dataframe to numpy array
                haplotypes = allel.HaplotypeArray(arr)   
                geno = haplotypes.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
                ac= geno.count_alleles()
                df1 = pd.read_pickle(y)#, index_col=0)
                df1 = df1.loc[(df1['POS'] >=  mini ) & (df1['POS'] <= maxi)]
                new_dataframe1 = df1.filter(['POS'])  # extract the position from the dataframe
                position1 = pd.DataFrame(new_dataframe1, columns=['POS']).to_numpy()
                position1 = position1.flatten()    # alter position to one dimensional numpy array for correct input for sequence diversity function
                df1 = df1.drop(['POS'], axis=1)
                arr1 = pd.DataFrame(df1).to_numpy()     # haplotype dataframe to numpy array
                haplotypes1 = allel.HaplotypeArray(arr1)   
                geno1 = haplotypes1.to_genotypes(ploidy=2)   # reshape haplotype array to genotype array for subsequent correct function input
                ac1= geno1.count_alleles()
                num,dem= allel.hudson_fst(ac,ac1) ## compute fst using allel package 
                fst = np.sum(num) / np.sum(dem)
                fst_lst=fst.tolist() ## tranform output array to a list in order to append it to the output list 
                output.append(fst_lst)
                df_fst = pd.DataFrame (output) ## trasnfom the output list to a data frame for easier reading  
                #fst_df_final=df_fst
        ## transpose the data frame 
        df_fst.columns = newData
        a=df_fst.T
        a.columns=['Fixation index']
        #df2 = fst_df_final.rename(columns=lambda x: x.strip('.pkl'))## name the data frame columns with the combinations calculated earlier 
        return a
