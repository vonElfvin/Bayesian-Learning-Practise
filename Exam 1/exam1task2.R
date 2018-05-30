#### Exam 1 Task 2: Modeling count data

### Help Functions

GibbsMixPois <- function(x, nComp, alpha, alphaGamma, betaGamma, xGrid, nIter){
  
  # Gibbs sampling for a mixture of Poissons
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   x - vector with data observations (counts)
  #   nComp - Number of mixture components to be fitted
  #   alpha - The prior on the mixture component weights is w ~ Dir(alpha, alpha,..., alpha) 
  #   alphaGamma and betaGamma - 
  #              The prior on the mean (theta) of the Poisson mixture components is 
  #              theta ~ Gamma(alphaGamma, betaGamma) [rate parametrization of the Gamma dist]
  #   xGrid - the grid of data values over which the mixture is evaluated and plotted
  #   nIter - Number of Gibbs iterations
  #
  # OUTPUTS:
  #   results$wSample     - Gibbs sample of mixture component weights. nIter-by-nComp matrix
  #   results$thetaSample - Gibbs sample of mixture component means.   nIter-by-nComp matrix
  #   results$mixDensMean - Posterior mean of the estimated mixture density over xGrid.
  
  
  ####### Defining a function that simulates from a Dirichlet distribution
  rDirichlet <- function(param){
    nCat <- length(param)
    thetaDraws <- matrix(NA,nCat,1)
    for (j in 1:nCat){
      thetaDraws[j] <- rgamma(1,param[j],1)
    }
    thetaDraws = thetaDraws/sum(thetaDraws) # Diving every column of ThetaDraws by the sum of the elements in that column.
    return(thetaDraws)
  }
  
  # Simple function that converts between two different representations of the mixture allocation
  S2alloc <- function(S){
    n <- dim(S)[1]
    alloc <- rep(0,n)
    for (i in 1:n){
      alloc[i] <- which(S[i,] == 1)
    }
    return(alloc)
  }
  
  # Initial values for the Gibbs sampling
  nObs <- length(x)
  S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
  theta <- rep(mean(x), nComp) # Each component is initialized at the mean of the data
  
  # Setting up the grid where the mixture density is evaluated.
  mixDensMean <- rep(0,length(xGrid))
  effIterCount <- 0
  
  # Setting up matrices to store the draws
  wSample <- matrix(0, nIter, nComp)
  thetaSample <- matrix(0, nIter, nComp)
  probObsInComp <- rep(NA, nComp)
  
  # Setting up the priors - the same prior for all components
  alpha <- rep(alpha, nComp) 
  alphaGamma <- rep(alphaGamma, nComp) 
  betaGamma <- rep(betaGamma, nComp) 
  
  # HERE STARTS THE ACTUAL GIBBS SAMPLING
  
  for (k in 1:nIter){
    message(paste('Iteration number:',k))
    alloc <- S2alloc(S) # Function that converts between different representations of the group allocations
    nAlloc <- colSums(S)
    
    # Step 1 - Update components probabilities
    w <- rDirichlet(alpha + nAlloc)
    wSample[k,] <- w
    
    # Step 2 - Update theta's in Poisson components
    for (j in 1:nComp){
      theta[j] <- rgamma(1, shape = alphaGamma + sum(x[alloc == j]), rate = betaGamma + nAlloc[j])
    }
    thetaSample[k,] <- theta
    
    # Step 3 - Update allocation
    for (i in 1:nObs){
      for (j in 1:nComp){
        probObsInComp[j] <- w[j]*dpois(x[i], lambda = theta[j])
      }
      S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
    }
    
    # Computing the mixture density at the current parameters, and averaging that over draws.
    effIterCount <- effIterCount + 1
    mixDens <- rep(0,length(xGrid))
    for (j in 1:nComp){
      compDens <- dpois(xGrid, lambda = theta[j])
      mixDens <- mixDens + w[j]*compDens
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
  }
  return(results = list(wSample = wSample, thetaSample = thetaSample, mixDensMean = mixDensMean))
}

### Setup
load(file = 'bids.Rdata')
bids.count = table(bids)

### Implementation
# a)
# prior
alpha0 = 1
beta0 = 1

# gamma is a conjugate prior to poisson, gamma -> Gamma(alpha0 + sum(bids), beta0 + n)
n = length(bids)
theta.grid = seq(3.3, 4, by=0.001)
theta.post = dgamma(theta.grid, shape = alpha0 + sum(bids), rate = beta0 + n)
plot(theta.grid, theta.post)

# b)
# calculate data density/distribution
data.distribution = bids.count / sum(bids.count)

# calculate the poisson density/distribtuion using the posterior mean
x.grid = seq(0,12, by=1)
# pois.distribution = dpois(x.grid, lambda=theta.grid[which.max(theta.post)])
pois.distribution = dpois(x.grid, lambda=mean(bids))

# plot and compare the two densities
plot(data.distribution, type="l", ylim=c(0,0.25))
lines(x.grid, pois.distribution)
# Comment: awful fit

# c)
alpha = 1 # uniform prior for weights
alphaGamma = 1
betaGamma = 1
nIterations = 500
x.grid = seq(0, 12)

gibbsMixPois.res.k2 = GibbsMixPois(bids, nComp=2, alpha=1, alphaGamma=1, betaGamma=1, xGrid=x.grid, nIter=nIterations)
gibbsMixPois.res.k3 = GibbsMixPois(bids, nComp=3, alpha=1, alphaGamma=1, betaGamma=1, xGrid=x.grid, nIter=nIterations)

# d)
lines(gibbsMixPois.res.k2$mixDensMean)
lines(gibbsMixPois.res.k3$mixDensMean)
# The two mixture distributions give similar fit and much closer to the data distribution
# I would choose k=2 since its simpler
