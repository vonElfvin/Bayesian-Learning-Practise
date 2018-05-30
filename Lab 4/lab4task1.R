#### Lab 4 Task 1: Poisson regression - the MCMC way.
### Libraries
library(mvtnorm)

### Functions
logPostPoisson = function(beta, X, y, mu.0, Sigma.0) {
  
  # likelihood
  log.likelihood = sum(y * X %*% beta - exp(X %*% beta))
  if(abs(log.likelihood)==Inf) log.likelihood = -20000
  
  # prior
  log.prior = dmvnorm(beta, mean = mu.0, sigma = Sigma.0, log = TRUE)
  
  # posterior
  log.posterior = log.likelihood + log.prior
  
  return(log.posterior)
}

Metropolis = function(theta, c, propSigma, logPostFunc, nBurnin, nSamples, ...) {
  
  # setup sampling
  p = length(theta)
  theta.samples = matrix(NA, nSamples, p)
  
  # initial values 
  theta.c = theta
  
  for(i in -nBurnin:nSamples) {
    # propose sample
    theta.p = as.vector(rmvnorm(1, mean = theta.c, sigma = c * propSigma))
    
    # calculate posterior likelihoods
    log.post.p = logPostFunc(theta.p, ...)
    log.post.c = logPostFunc(theta.c, ...)
    
    # calculate alpha
    alpha = min(1, exp(log.post.p - log.post.c))
    
    # with probabiliy alpha, choose proposed sample as sample
    u = runif(1, min = 0, max = 1)
    if(u < alpha) theta.c = theta.p
    
    # save sample if not burnin phase
    if(i>0) theta.samples[i,] = theta.c
  }
  
  return(theta.samples)
}

### Setup
data.ebay = read.table("eBayNumberOfBidderData.dat", header=TRUE)
X = as.matrix(data.ebay[,-1])
y = data.ebay[,1]
n = length(y)
p = dim(X)[2]

### Implementation
# a)
glm.fit = glm(nBids ~ 0 + ., data=data.ebay, family="poisson")
summary(glm.fit)
#   Coefficients:
#               Estimate Std. Error z value  Pr(>|z|)    
#   Const        1.07244    0.03077  34.848   < 2e-16 *** #significant
#   PowerSeller -0.02054    0.03678  -0.558    0.5765    
#   VerifyID    -0.39452    0.09243  -4.268 0.0000197 *** #significant
#   Sealed       0.44384    0.05056   8.778   < 2e-16 *** #significant
#   Minblem     -0.05220    0.06020  -0.867    0.3859    
#   MajBlem     -0.22087    0.09144  -2.416    0.0157 *   #significant
#   LargNeg      0.07067    0.05633   1.255    0.2096    
#   LogBook     -0.12068    0.02896  -4.166 0.0000309 *** #significant
#   MinBidShare -1.89410    0.07124 -26.588   < 2e-16 *** #significant
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# b)
# prior
mu.0 = rep(0, p)
Sigma.0 = 100 * solve(t(X) %*% X)

# optim
optim.res = optim(mu.0, logPostPoisson, X, y, mu.0, Sigma.0, 
                  gr = NULL, method = "BFGS", control = list(fnscale = -1), hessian = T)

# extract results
beta.mode = optim.res$par
beta.neg.inv.hessian = -solve(optim.res$hessian)

# c)
nBurnin = 1000
nSamples = 4000
c = 0.5
beta.samples = Metropolis(mu.0, c, beta.neg.inv.hessian, 
                            logPostPoisson, nBurnin, nSamples,  
                            X, y, mu.0, Sigma.0)

# d)
# example auction
auction = c(1, 1, 1, 1, 0, 0, 0, 1, 0.5)

# sample nBids
nBids.samples = numeric()
for(i in 1:nSamples) {
  nBids.samples[i] = rpois(1, exp(auction %*% beta.samples[i,]))
}
# plot samples
barplot(table(nBids.samples))

# calculate probabilty of no bidders
prob = sum(ifelse(nBids.samples == 0, 1, 0)) / nSamples # 0.3825
