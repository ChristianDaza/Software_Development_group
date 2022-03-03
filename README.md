# Genome Analytics
![](https://img.shields.io/badge/licence-Genome%20Analytics-blue) ![](https://img.shields.io/badge/version-1.0.0-blue) ![](https://img.shields.io/badge/platforms%20-macOS--64%20%7C%20win--64-lightgrey)

##  SNP Browser  Software

SNP Browser is a fast and accurate software that allows for genomic data retrieval of SNPs located on chromosome 22 in conjunction with genotype and allele frequencies across populations. A variety of population summary statistics, which include, Nucleotide diversity, Homozygosity, Tajima’s D, and fixation index (FST) can be computed, which are returned in a tabular format and plotted along sliding windows, with the additional opportunity to download the results as a text file. 


SNP Browser was developed by Genome Analytics, a group composed of four MSc Bioinformatics students from the Queen Mary University of London under the guidance of Professor Conrad Bessant and Dr Matteo Fumagalli.


- Documentation:  https://github.com/ChristianDaza/Software_Development_group/blob/Main/Documentation
- Source code: https://github.com/ChristianDaza/Software_Development_group/tree/Main/website 


## Getting  started



### Running our Software
Our software can be run in the command line, Visual studio Code and PyCharm, but we will follow the command line approach.


### Installing
To run our software from your local machine, use the command line to navegate into your desktop directory and download this repository using the following code:

```
git clone https://github.com/ChristianDaza/Software_Development_group.git
```

### Required packages
Using the command line install the following packages, found in the [requirement](https://github.com/ChristianDaza/Software_Development_group/blob/Main/requirements.txt) file.

- Python 3.9.5
- Flask==2.0.2
- Flask-WTF==1.0.0
- WTForms==3.0.1
- Flask_Bootstrap==3.3.7.1
- Flask_APScheduler==1.12.3
- Pandas==1.4.0
- numpy==1.20.3
- scikit-allel==1.3.5
- plotly==5.5.0

Use the following command to install all the required packages.

```
pip install -r requirements.txt 
```

### Runing SNP Browser
Within the command line  move into the  website directory and type the following command. 

``` 
python3 application.py 
```

Then copy and paste the URL into your browser (we recommend Google Chrome or Safari), which will direct you to the homepage of SNP Browser, where you can start your groundbreaking research.. 

## Authors 
- Christian David Arenas Daza: [ChristianDaza](https://github.com/ChristianDaza)
- Eri Krasniqi: [erikras12](https://github.com/erikras12)
- Rym Khadidja Lebboukh: [RYM213](https://github.com/RYMY213)                                
- Yasemin Bridges: [yaseminbridges](https://github.com/yaseminbridges)   

## Development status
Please note that SNP Browser is no longer under the development phase and will only be maintained.

## How to cite our software

APA:
```
Arenas Daza, C. D., Krasniqi, E., Lebboukh, R. K., & Bridges, Y. (2022). SNP Browser (Version 1.0.0) [Computer software]. https://github.com/ChristianDaza/Software_Development_group
```

BibTeX:
```
@software{Arenas_Daza_SNP_Browser_2022,
author = {Arenas Daza, Christian David and Krasniqi, Eri and Lebboukh, Rym Khadidja and Bridges, Yasemin},
doi = {10.5281/zenodo.1234},
month = {3},
title = {{SNP Browser}},
url = {https://github.com/ChristianDaza/Software_Development_group},
version = {1.0.0},
year = {2022}
}
```

## Acknowledgments
-  We would like to thank Professor Conrad Bessant and Dr Matteo Fumagalli for the support and advice provided during the development of this software.


