#Before running this script, I went to https://sourceforge.net/projects/mcmc-jags/files/latest/download
# and downloaded JAGS sampler. 

#Then I looked for MixSIAR in the cran library.
install.packages("MixSIAR")

#After the package is installed, I need to load it into my computer memory
library(MixSIAR)
#load the library for rjags
library(rjags)
library(R2jags)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(reshape)
library(reshape2)
library(lattice)
library(MCMCpack)
library(ggmcmc)
library(coda)
library(loo)
library(bayesplot)

install.packages("splancs")
library(splancs)
update.packages(ask = FALSE, checkBuilt = TRUE)




