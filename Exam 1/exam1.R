##### Exam 1

#### Exam 1 Task 1: Bayesian inference for the rice distribution

### Help Functions

# Random number generator for the Rice distribution
rRice <-function(n = 1, theta = 1, psi = 1){
  x <- rnorm(n = n, mean = 0, sd = sqrt(psi))
  y <- rnorm(n = n, mean = theta, sd = sqrt(psi))
  return(sqrt(x^2+y^2))
}

# log density of rice distribution
logdrice = function(x, psi, theta) {
  return(log(x/psi) - (x^2 + theta^2) / (2 * psi) + log(besselI(x * theta / psi, nu=0)))
}

### Setup
riceData <- c(1.556, 1.861, 3.135, 1.311, 1.877, 0.622, 3.219, 0.768, 2.358, 2.056)
# psi = 1

# a)
# log posterior distribution of theta given rice model and a prior proportionate to c
logPostRice = function(theta, x) {
  # likelihood
  log.likelihood = sum(logdrice(x, psi=1, theta))
  
  # prior
  log.prior = 0 # assumed to be proprotionate to constant
  
  # posterior
  log.posterior = log.likelihood + log.prior
  
  return(log.posterior)
}

grid.width = 0.01
theta.grid = seq(0, 3, by=grid.width)
logPost = numeric()

for(i in 1:length(theta.grid)) {
  logPost[i] = logPostRice(theta.grid[i], riceData)
}

plot(theta.grid, exp(logPost) / (sum(exp(logPost)) * grid.width), type="l", lwd=2)

# b)
# optim
optim.res = optim(mean(riceData), logPostRice, riceData,
                  gr = NULL, method = "L-BFGS-B", lower = 0.001,
                  control = list(fnscale=-1), hessian = TRUE)
# extract mode and variance
theta.mode = optim.res$par
theta.inv.neg.hessian = as.numeric(-solve(optim.res$hessian))

# line the result
lines(theta.grid, dnorm(theta.grid, mean=theta.mode, sd=sqrt(theta.inv.neg.hessian)), col="red", lwd=2)

# c)
nSamples = 10000
xtilde.samples = numeric()
theta.samples = rnorm(nSamples, mean=theta.mode, sd=sqrt(theta.inv.neg.hessian))
for(i in 1:nSamples) {
  xtilde.samples[i] = rRice(n = 1, theta = theta.samples[i], psi = 1)
}
# plot(density(xtilde.samples))
hist(xtilde.samples, 100, freq=FALSE)

#### Exam 1 Task 2: Modeling count data

### Help Functions

# Code for Problem 3 - Exam in Bayesian Learning 2017-05-30
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
load(file = 'bids.RData')    # Loading the vector 'bids' into workspace
bids.count = table(bids)

# a)
grid.width = 0.001
theta.grid = seq(3.4, 3.9, by=grid.width)
alpha.prior = 1
beta.prior = 1
n = length(bids)
x.mean = mean(bids)
plot(theta.grid, dgamma(theta.grid, alpha.prior + n * x.mean, beta.prior + n), type="l", lwd=2,
     xlab="theta", ylab="density")

# b)
x.grid = seq(min(bids), max(bids))
nSamples = 10000
theta.samples = rgamma(nSamples, alpha.prior + n * x.mean, beta.prior + n)
poisson.densities = matrix(NA, nSamples, length(x.grid))

for(i in 1:nSamples) {
  poisson.densities[i,] = dpois(x.grid, lambda = theta.samples[i])
}
mean.pois.dens = apply(poisson.densities, 2, mean)
data.dens = bids.count / sum(bids.count)
plot(x.grid, mean.pois.dens, lwd=2, col="red", type="l")
lines(x.grid, data.dens, ylim=c(0, 0.3), lwd=2)

# c)
gibbs.res.k2 = GibbsMixPois(bids, nComp=2, alpha = 1, alphaGamma = 1, betaGamma = 1, x.grid, nIter=500)
gibbs.res.k3 = GibbsMixPois(bids, nComp=3, alpha = 1, alphaGamma = 1, betaGamma = 1, x.grid, nIter=500)

# d)
plot(x.grid, data.dens, ylim=c(0, 0.3), lwd=2, type="l")
lines(x.grid, gibbs.res.k2$mixDensMean, col="blue", lwd=2)
lines(x.grid, gibbs.res.k3$mixDensMean, col="green", lwd=2)

#### Exam 1 Task 3: Regression
### Help Functions
# Defining a function that simulates from the scaled inverse Chi-square distribution
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){
  # Direct sampling from a Gaussian linear regression with conjugate prior:
  #
  # beta | sigma2 ~ N(mu_0, sigma2*inv(Omega_0))
  # sigma2 ~ Inv-Chi2(v_0,sigma2_0)
  # 
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   y - n-by-1 vector with response data observations
  #   X - n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
  #   mu_0 - prior mean for beta
  #   Omega_0  - prior precision matrix for beta
  #   v_0      - degrees of freedom in the prior for sigma2
  #   sigma2_0 - location ("best guess") in the prior for sigma2
  #   nIter - Number of samples from the posterior (iterations)
  #
  # OUTPUTS:
  #   results$betaSample     - Posterior sample of beta.     nIter-by-nCovs matrix
  #   results$sigma2Sample   - Posterior sample of sigma2.   nIter-by-1 vector
  
  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  betaHat <- solve(XX,t(X)%*%y)
  Omega_n = XX + Omega_0
  mu_n = solve(Omega_n,XX%*%betaHat+Omega_0%*%mu_0)
  v_n = v_0 + n
  sigma2_n = as.numeric((v_0*sigma2_0 + ( t(y)%*%y + t(mu_0)%*%Omega_0%*%mu_0 - t(mu_n)%*%Omega_n%*%mu_n))/v_n)
  invOmega_n = solve(Omega_n)
  
  # The actual sampling
  sigma2Sample = rep(NA, nIter)
  betaSample = matrix(NA, nIter, nCovs)
  for (i in 1:nIter){
    
    # Simulate from p(sigma2 | y, X)
    sigma2 = rScaledInvChi2(n=1, df = v_n, scale = sigma2_n)
    sigma2Sample[i] = sigma2
    
    # Simulate from p(beta | sigma2, y, X)
    beta_ = rmvnorm(n=1, mean = mu_n, sigma = sigma2*invOmega_n)
    betaSample[i,] = beta_
    
  }
  return(results = list(sigma2Sample = sigma2Sample, betaSample=betaSample))
}

### Setup
# library
library(mvtnorm)

# data
load("cars.RData")
y = cars[,1]
X = as.matrix(cars[,-1])
cov.names = names(cars[,-1])

# priors
mu.0 = rep(0,4)
Omega.0 = 0.01 * diag(4)
v.0 = 1
sigma2.0 = 36

# a)
nIter = 1000
blr.res = BayesLinReg(y = y, 
                      X = X, 
                      mu_0 = mu.0, 
                      Omega_0 = Omega.0, 
                      v_0 = v.0,
                      sigma2_0 = sigma2.0, 
                      nIter = nIter)

beta.samples = blr.res$betaSample
colnames(beta.samples) = cov.names
sigma2.samples = blr.res$sigma2Sample

# plot marginal distributions
hist(beta.samples[,"intercept"], 30, freq=FALSE)
hist(beta.samples[,"weight"], 30, freq=FALSE)
hist(beta.samples[,"sixcyl"], 30, freq=FALSE)
hist(beta.samples[,"eightcyl"], 30, freq=FALSE)

# point estimates: linear loss -> posterior mean
beta.estimates = apply(beta.samples, 2, median)
# intercept    weight    sixcyl  eightcyl 
# 33.780135 -3.117405 -4.262106 -6.151775 

# calculate equal tail probability intervals
beta.intervals = apply(beta.samples, 2, function(x) quantile(x, prob=c(0.025, 0.975)))
#       intercept    weight    sixcyl  eightcyl
# 2.5%   29.80489 -4.634012 -7.179988 -9.618449
# 97.5%  37.57732 -1.612022 -1.268775 -2.632720
# no 0 included in any interval

# b)
beta.diff = beta.samples[,"sixcyl"] - beta.samples[,"eightcyl"]
CI.diff = quantile(beta.diff, prob=c(0.025, 0.975))
#       2.5%      97.5% 
# -0.8533905  4.7568483 
# 0 included -> no difference
hist(beta.diff, freq=FALSE) # sixcyl tend to have higher mpg, however not guaranteed

# c)
car = c(1, 3.5, 0, 0)

# sampling
mpg.samples = numeric()
nSamples = 1000
for(i in 1:nSamples) {
  mpg.samples[i] = car %*% beta.samples[i,] + rnorm(1, mean = 0, sd = sqrt(sigma2.samples[i]))
}

# plot density
hist(mpg.samples, 30, freq=FALSE)

#### Exam 1 Task 4: Geometric data and decisions