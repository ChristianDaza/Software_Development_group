# Genome Analytics

##  SNP Browser  Software

SNP browser is a fast and accurate application that allows information retrieval of SNPs from Chromosome 22 in conjunction with genotype and allele frequency.This application also allows the user to  perform genetic statistical testing; Nucleotide diversity, homozygosity, Tajima's D, and fixation index. The return results are display, plotted, and can be downloaded as a txt file. 

This software was developed by the Genomic Analytics group (4 MSc Bioinformatics students From Queen Mary University) under the guidance of Professor Conrad Bessant and Dr Matteo Fumagalli, in response to the need for an inclusive platform that not only allows the querying and extraction of genomic information but also statistical testing.

- Documentation:  https://github.com/ChristianDaza/Software_Development_group/blob/Main/website/Documentation
- Source code: https://github.com/ChristianDaza/Software_Development_group/tree/Main/website
- Meet our team: 



## Getting  started - End user

### The following packages will be required  for this application 


- Python 3.9.5
- Flask==2.0.2
- Werkzeug==2.0.3
- Flask-WTF==1.0.0
- WTForms==3.0.1
- Flask_Bootstrap==3.3.7.1
- Flask_APScheduler==1.12.3
- reuqest==2.27.1
- Pandas==1.4.0
- scikit-allel==1.3.5


To install the above packages:
```
pip install Flask-WTF  
pip install Werkzeug
pip install WTForms
pip install Flask-Bootstrap
pip install Flask-APScheduler
pip install requests
pip install scikit-allel
```

### Installing
To run this website from your local machine please download the "website" directory

### Runing our application 

Our application can be run various ways including code editors such as Visual Studio Code, Pycharm or in the command line.
In all cases open or navegate into the website directory in your local machine, and ran the following command.

``` python3 application.py ```

The either click or copy and paste the URL into your browser, we recommend Google Chrome. 

## Citation and reference

### How shoudl be cited in a document 

#### Licence

#### Contributing







The aim of this project is to create and deployed a piece of software with the following specifications :

The user should be able to retrieve SNP information given either a genomic coordinate (chromosome, start and end), SNP name (rs value), or gene name (or any aliases associated to it).

The application should return the following information for each SNP: name (rs value), genomic position, genotype frequencies, and allele frequency. Frequencies should be provided for each population separately.

If multiple SNPs are returned, the user should be able to select the population(s) and summary statistics of interest, and the application will calculate them and plot their distribution in slidingwindows along the region. The user will also be able to download a text file with the values of summary statistics. At least three summary statistics should be reported, for instance: one index of genetic diversity, one on haplotype diversity (e.g. homozygosity), one test against neutrality (e.g. Tajima’s D). Additionally, if multiple populations are selected, then population genetic variation (FST value) for each pair of populations should be reported.
