import allel
import pandas as pd


####### gene names retrieved from annotated vcf file #######

# read in annotated SNP file
callset = allel.read_vcf('annotated_chr22.vcf', fields='ANN', transformers=allel.ANNTransformer())
# extract only the gene names
genes = callset['variants/ANN_Gene_Name'] 
# turn into data frame
df = pd.DataFrame(data=genes, columns = ['Genes'])

####### read in NCBI gene data and map only dataframe to retrieve corresponding gene names #######

NCBI_genes = pd.read_csv('Homo_sapiens.gene_info', delimiter='\t')
# isloate chromosome 22 genes
NCBI_genes = NCBI_genes[NCBI_genes.chromosome=='22']
# retrieve only gene symbol and synonyms
NCBI_genes_aliases = NCBI_genes[['Symbol','Synonyms']]
# create dictionary for mapping
genes_aliases_dict = dict(zip(NCBI_genes_aliases.Symbol, NCBI_genes_aliases.Synonyms))
# map aliases
dfseries=pd.Series(df['Genes'])
df2 = dfseries.map(genes_aliases_dict)
# create final dataframe
df['Aliases']=df2
# to csv
df.to_csv('GENES_AND_ALIASES_FINAL.csv')








