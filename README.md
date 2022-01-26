# Software_Development_group

The aim of this project is to create and deployed a piece of software with the following specifications :

The user should be able to retrieve SNP information given either a genomic coordinate (chromosome, start and end), SNP name (rs value), or gene name (or any aliases associated to it).

The application should return the following information for each SNP: name (rs value), genomic position, genotype frequencies, and allele frequency. Frequencies should be provided for each population separately.

If multiple SNPs are returned, the user should be able to select the population(s) and summary statistics of interest, and the application will calculate them and plot their distribution in slidingwindows along the region. The user will also be able to download a text file with the values of summary statistics. At least three summary statistics should be reported, for instance: one index of genetic diversity, one on haplotype diversity (e.g. homozygosity), one test against neutrality (e.g. Tajima’s D). Additionally, if multiple populations are selected, then population genetic variation (FST value) for each pair of populations should be reported.
