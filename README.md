# Nonparametric Double Robustness R Code

This in progress R code in this repository is being used to update results for the "Nonparametric Double Robustness Paper," currently in the arXiv: https://arxiv.org/abs/1711.07137. It was run from the command line using the following:

Rscript --no-save --no-restore --verbose npDR_version2.R 1 10000 > outV2.Rout 2>&1

where the "1" and "10000" variables represent the seed and the number of Monte Carlo simulations, respectively.

The npDR_version2.R file runs the simulations, and outputs a results.txt file. The post.R takes the results and generates figures and tables found in the paper. 


