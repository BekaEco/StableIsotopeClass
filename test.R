library(MixSIAR)
library(rjags)
library(R2jags)
library(splancs)


source("MixSIAR_readDF_functions.r")

set.seed(850) #make any given result repeatable...

# Specify means and variances for sources and mixtures ###############################
# establish source means (4 sources (rows), delC in 1st column, delN in 2nd column5
src_means <- matrix(c(-20,5,-20,15,-10,15,-10,5), nrow = 4, ncol=2, byrow=TRUE)

# establish source variances (4 sources (rows), delC in 1st column, delN in 2nd column)
src_sds <- matrix(c(0.6,0.6,0.6,0.6,1.12,1.12,1.12,1.12), nrow = 4, ncol=2, byrow=TRUE)

#establsh mixture means (delC in 1st column, delN in 2nd column)
mix_means <- c(-15,10)

#establsh mixture variances (delC in 1st column, delN in 2nd column)
mix_sds <- c(0.6,1.12)

# Simulate data #####################################################################
N <- 5 #sample size
num.src <- 4 #number of sources
num.iso <- 2 # number of tracers
S <- array(0, dim=c(num.src,num.iso,N))
X <- array(0, dim=c(num.iso,N))

for(i in 1:N) {
  for(src in 1:num.src) {
    for(iso in 1:num.iso) {
      S[src,iso,i] <- rnorm(1, mean = src_means[src,iso], sd = src_sds[src,iso] ) #source data
      X[iso,i] <- rnorm(1, mean = mix_means[iso], sd = mix_sds[iso] )		    #mixture data
    }
  }
}

#PLOT THE DATA #######################################################################
plot(S[1,1,],S[1,2,],ylim=c(1,19),xlim=c(-23,-7),pch =0, xlab="Isotope 1 (per mil)", ylab="Isotope 2 (per mil)") 	#first source
points(S[2,1,],S[2,2,],pch=2) 								#second source
points(S[3,1,],S[3,2,],pch=5) 								#third source
points(S[4,1,],S[4,2,],pch=6) 								#fourth source
points(X[1,],X[2,],pch = 16)							  	# mixture

#Convert list to dataframe (bogus)
src1 <- t(S[1,,])
src2 <- t(S[2,,])
src3 <- t(S[3,,])
src4 <- t(S[4,,])
sources <- data.frame(Source = rep(1:num.src, each=N), rbind(src1, src2, src3, src4))

mixtures <- data.frame(X1=X[1,], X2=X[2,])


#Run MixSIAR #######################################################################
mix <- load_mix_data(filename=mixtures,
                     iso_names=c("X1","X2"),
                     factors=NULL,
                     fac_random=NULL,
                     fac_nested=FALSE,
                     cont_effects=NULL)

## ------------------------------------------------------------------------
source <- load_source_data(filename=sources,
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="raw",
                           mix)

## ------------------------------------------------------------------------
disc <- data.frame(1:4, MeanX1=rep(0,num.src ), SDX1=rep(0,num.src ), MeanX2=rep(0,num.src ), SDX2=rep(0,num.src ))

## ------------------------------------------------------------------------
discr <- load_discr_data(disc, mix)

## -------------------------------------------------------------
  # Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

## ------------------------------------------------------------------------
  # Calculate the convex hull area, standardized by source variance
  polygonArea <- calc_area(source=source,mix=mix,discr=discr)

## -------------------------------------------------------------
  # default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
  plot_prior(alpha.prior=1,source)

## -------------------------------------------------------------
  # Write the JAGS model file
  model_filename <- "MixSIAR_model.txt"
  resid_err <- TRUE
  process_err <- FALSE
  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

## -------------------------------------------------------------
  jags.1 <- run_model(run="test", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)

## -------------------------------------------------------------
  jags.1 <- run_model(run="normal", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)

## -------------------------------------------------------------
  output_JAGS(jags.1, mix, source)
  
  
  #PLOT THE DATA #######################################################################
  plot(S[1,1,],S[1,2,],ylim=c(1,19),xlim=c(-23,-7),pch =0, xlab="Isotope 1 (per mil)", ylab="Isotope 2 (per mil)") 	#first source
  points(S[2,1,],S[2,2,],pch=2) 								#second source
  points(S[3,1,],S[3,2,],pch=5) 								#third source
  points(S[4,1,],S[4,2,],pch=6) 								#fourth source
  points(X[1,],X[2,],pch = 16)							  	# mixture  
