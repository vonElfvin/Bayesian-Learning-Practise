#### Exam 3 Task 2: Regression

### Libraries
library(MASS)
library(mvtnorm)

### Setup
BostonHousing = Boston
y = BostonHousing$medv
X = cbind(1,BostonHousing[,1:13]) # Adding a column of ones for the intercept
names(X)[1] <- "intercept"
covNames <- names(X)
y <- as.numeric(y)
X <- as.matrix(X)
p = dim(X)[2]

### Help Functions
# Defining a function that simulates from the scaled inverse Chi-square distribution
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

# a)
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

# perform sampling
nSamples = 5000
res = BayesLinReg(y, X, rep(0,p), 10^2*diag(p), 1, 6^2, nSamples)

# extract result
sigma2.samples = res$sigma2Sample
beta.samples = res$betaSample
colnames(beta.samples) = covNames

# extract beta of interest
beta.lstat = beta.samples[,"lstat"]

# calculate HPD: symmetric distirbution (student-T) -> HPD same as equal-tail
HPDI.lstat = quantile(beta.lstat, prob=c(0.025, 0.975))
#       2.5%      97.5% 
# -0.6363353 -0.4310156 

# b)
# setup house
house.before = X[9,]
house.after = house.before
house.after["lstat"] = house.after["lstat"]*0.7

# perform predictions sampling
prices.before = numeric()
prices.after = numeric()
for(i in 1:nSamples) {
  prices.before[i] = beta.samples[i,] %*% house.before + rnorm(1, mean=0, sd=sqrt(sigma2.samples[i]))
  prices.after[i] = beta.samples[i,] %*% house.after + rnorm(1, mean=0, sd=sqrt(sigma2.samples[i]))
}

# hist data
# prices.change = prices.after - prices.before
prices.change = beta.samples[,p]*(house.after["lstat"]-house.before["lstat"])

# HPD price change
HPDI.price.change = quantile(prices.change, prob=c(0.025, 0.975))
