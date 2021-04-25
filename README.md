
# TOP

################################################
# Targeted selection by integrating organismal and molecular traits through predictive analytics
################################################

# Project Introduction:
1. The likelihood of reduced yields of major food crops due to the changing climate continues to rise, as does the global population, thus making the development of genetically improved, climate-resilient crops a research priority. Genomic prediction and selection have been implemented in many crops to accelerate the breeding process in public and private breeding programs. Plant breeders are usually interested in improving multiple traits, but breeding for two or more traits simultaneously is generally more difficult than breeding for a single trait. However, genomic prediction can be advanced by shifting the focus from genomic-predicted values to identification of individual plants that come closest to an ideotype, or target variety, which combines merits from multiple traits. 
2. We present an integrative multi-trait breeding strategy that incorporates agronomic and omics traits (transcriptomes and metabolomes) to predict the best performing candidates to create new varieties through a machine learning algorithm. This algorithm, called Target-Oriented Prioritization (TOP), learns the inherent correlations among traits in a training population, balances the selection of multiple traits simultaneously, and predicts the degree of similarity between an untested genotype and a target, which can be a current commercial variety.
3. We examined this strategy in a maize NCII population with 5820 F1 hybrids derived from crossing 194 maternal CUBIC lines and 30 paternal lines, and calculated the accuracy of identifying a breeding candidate of a predefined target. This strategy was further extrapolated to two independent maize populations of diverse inbred lines and a rice population of recombinant inbred lines. 
4. The genotype (156,269 SNPs) and phenotypes (18 agronomic traits) of the 5,820 F1 hybrids were deposited in this branch of GitHub and in the MAIZEGO resource website (http://www.maizego.org/Resources.html). Genotype data consists of genotype indicators. The value 0 is for the homozygote of the major allele, the value 1 is for the heterozygote and the value 2 is for the homozygote of the minor allele. A total of 18 agronomic traits including flowering traits (days to tassel, days to anther and days to silk), plant architecture traits (plant height, ear height, ear leaf width, ear leaf length, tassel length, tassel branch number) and yield traits (cob weight, ear weight, ear diameter, ear length, ear row number, kernel number per ear, kernel number per row, kernel weight per ear, length of baren tip).
5. The R script of TOP algorithm and related demo data are provided.

# Genotypic and phenotypic data
We used four datasets for testing TOP algorithm, including maizeNCII, maize368, maize282 and rice210. The maizeNCII data are provided here in the GitHub, other three datasets are publicly available data which can be accessed from their original publications.

# TOP algorithm and demos
1. The "TOP_Demo.r" is a user-friendly R code which only needs to install rrBLUP package for genome prediction, a customed R script for building TOP model and making evaluation.
2. We provide three R scripts: TOP_Demo.r, GP.r and TOP.r.
All things need to do is run the "TOP_Demo.r", in which, the program will source the functions from "GP.r" and "TOP.r", 
and import two sample data "demo_geno.csv" and "demo_pheno.txt". The deom genotyic data contains 1619 bins of 210 individuals and the demo phenotypic data contains 4 traits of the 210 individuals. 
3. The function from GP.r can perfrom genomic predictions based on GBLUP model, which needed to be used as the inputs of TOP algorithm to learn the multi-trait inherent correlations.
Then, the function from TOP.r will build the TOP model by fitting the pairwised predicted and observed values for each individual, with parameterizing the multi-trait weights. The optimization of the TOP model will learn the optimal weights, which will be used to select the candidate individual from the testing pool, that shows the highest similarity on multiple traits. 









