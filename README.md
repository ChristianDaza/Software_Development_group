# Genome Analytics
![](https://img.shields.io/badge/licence-Genome%20Analytics-blue) ![](https://img.shields.io/badge/version-1.0.0-blue) ![](https://img.shields.io/badge/platforms%20-macOS--64%20%7C%20win--64-lightgrey)

##  SNP Browser  Software

SNP browser is a fast and accurate application that allows information retrieval of SNPs from Chromosome 22 in conjunction with genotype and allele frequency.This application also allows the user to  perform genetic statistical testing; Nucleotide diversity, homozygosity, Tajima's D, and fixation index. The return results are display, plotted, and can be downloaded as a txt file. 

This software was developed by the Genomic Analytics group (4 MSc Bioinformatics students From Queen Mary University) under the guidance of Professor Conrad Bessant and Dr Matteo Fumagalli, in response to the need for an inclusive platform that not only allows the querying and extraction of genomic information but also statistical testing.

- Documentation:  https://github.com/ChristianDaza/Software_Development_group/blob/Main/Documentation
- Source code: https://github.com/ChristianDaza/Software_Development_group/tree/Main/website 



## Getting  started

### The following packages will be required  for this application 


- Python 3.9.5
- Flask==2.0.2
- Werkzeug==2.0.3
- Flask-WTF==1.0.0
- WTForms==3.0.1
- Flask_Bootstrap==3.3.7.1
- Flask_APScheduler==1.12.3
- requests==2.27.1
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
pip install plotly
```

### Installing
To our software from your local machine please download the this repository using the following code:
```
git clone https://github.com/ChristianDaza/Software_Development_group.git
```

### Runing our application 

Our software can be run in the command line, Visual studio Code and PyCharm.

Using the command line  move into the  website directory and then type the following command 

``` python3 application.py ```

Then copy and paste the URL into your browser (we recommend Google Chrome or Safary), which then will direct you into the homepage of SNP Browser, where you can start your ground breaking research. 

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
-  We like to thank Professor Conrad Bessant and Dr Matteo Fumagalli for their support and advice provided during the developemnt of this software.


