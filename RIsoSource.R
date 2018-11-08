##############################################
#ISOSOURCE CODE
##############################################

############# USER INPUT #####################
#Necessary inputs:
# 1. A vector with source isotope values
# 2. A vector with mixture isotope values
# 3. A tolerance vector to establish the tolerance for each isotope (in per mil) for with the modeled mixture is a match to the observed.
# 4. An increment value (in proportion) for conducting trials across sources
# 5. A tolerance level to establish the tolerance for added proportions of all isotopes (e.g. 0.98-1.02)

# Example 1: Beka's mixture of blue, pink, and purple paint!
sources <- as.matrix(data.frame("SourceIso1" = c(4,5,6)))
consumer <- c(5)
toleranceVec <- c(0.1)
increment <- 0.01
tolerance.prop <- 0.01

#Example with two isotopes, four sources!
SourceIso1 <- c(-30,-25,-20, -30)
SourceIso2 <- c(10,10,5, 0)
sources <- cbind(SourceIso1,SourceIso2)
toleranceVec <- c(0.1, 0.2)
increment <- 0.02
tolerance.prop <- 0.01

##############################################
##############################################

############# FUNCTIONS ######################
all.options <- function(npars, min, max, increment, tolerance) {
    level.list <- list()

    for (i in 1:npars) {
        level.list[[i]] <- seq(from=min, to=max, by=increment)
    }
    grid <- expand.grid(level.list)
    sums <- rowSums(grid[,1:npars])
    grid <- grid[sums < (1 + tolerance) & sums > (1-tolerance), ]
    return(grid)
}

sumproduct <- function(x, sources){
    return(sum(x %*% sources))
}

isoSource <- function(sources, consumer, min=0, max=1, tolerance.prop, tolerance.iso, increment){
    
    #Determine number of sources and isotopes for the sources matrix
    numSources <- length(sources[,1])
    numIsos <- length(sources[1,])
    
    #Calculate matrix of all possible combinations of sources
    x <- all.options(npar=numSources, min=min,max=max,increment=increment, tolerance=tolerance.prop)
    numTrials <- length(x[,1])
    
    # Now create mixtures based on possible source proportions and isotopic ratios of the sources.
    mixture <- as.data.frame(matrix(data=NA, nrow=numTrials, ncol=numIsos), flag=F)
    for (i in 1:numIsos)  mixture[,i] <- apply(x,1,sumproduct, sources=sources[,i])
    
    # Now find the proportional contributions that solve the mass balance within the specified tolerance.
    test <- matrix(data=F, nrow=numTrials, ncol = numIsos)
    for (jj in 1:numIsos) test[,jj] <- abs(mixture[,jj] - consumer[jj]) <= tolerance.iso[jj]
    mixture$flag <- apply(test,1, all)
    
    #Final matrix of proptrional contributions from each source that satisfies the mass balance.
    validMixtures <- subset(data.frame(x,mixture), subset=mixture$flag==T)
    
    return(validMixtures)
}


##############################################
##############################################

############# MAIN ######################
#Run IsoSource code!
validMixture <- isoSource(sources=sources, consumer=consumer, min=0, max=1, tolerance.prop=tolerance.prop, tolerance.iso = toleranceVec, increment=increment)

#Plot stuff now!
long_DF <- validMixture %>% gather(Source, Proportion, Var1:Var3) # Need to change range of variables (VAR1:VAR3) for # sources
ggplot(long_DF, aes(x = Proportion)) + geom_histogram(binwidth=0.05)+facet_grid(~Source)+theme_bw() + scale_x_continuous(limits = c(0, 1))

#rm(mixture, test, x, SourceIso1, SourceIso2) #Clean house




