## readme.txt

1. System requirements

The "TOP_Demo.r" is a user-friendly R code which only needs to install rrBLUP package for genome prediction.

2. Intallation guide

Our code needs to install rrBLUP package, two R scripts containing customed functions can be easily and fastly loaded into R environment.

3 and 4. Demo and Intructions for use

We provide three R scripts: TOP_Demo.r, GP.r and TOP.r.
All things need to do is run the "TOP_Demo.r", in which, the program will source the functions from "GP.r" and "TOP.r", 
and import two sample data "RIL-Genotype.csv" and "data_trait.txt". 
The first data contains 1619 bins of 210 individuals and the second data contains 4 traits of the 210 individuals. 
The function from GP.r can obtain the trait predictions of training set and testing set using GBLUP.
Then, the function from TOP.r will construct TOP model, test TOP model, obtain the optimal weights and select similar individuals for a given target.









