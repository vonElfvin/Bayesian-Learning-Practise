#### Exam 2 Task 1: Bayesian inference for Cauchy data
library(mvtnorm)

### Help Functions

dCauchy <- function(x, theta = 0, gamma = 1){
  return(dens = (1/(pi*gamma))*(1/(1+((x-theta)/gamma)^2)))
}

dlognormal <- function(x, mu, sigma2){
  return(dens = (1/(sqrt(2*pi*sigma2)*x))*exp((-1/(2*sigma2))*(log(x)-mu)^2))
}

logPostCauchyNormalPrior = function(theta, gamma, y, mu0, sigma0) {
  # likelihood
  log.likelihood = sum(log(dCauchy(y, theta, gamma)))
  
  # prior
  log.prior = dnorm(theta, mu0, sd=sigma0, log = TRUE)
  
  # posterior
  log.posterior = log.likelihood + log.prior
  
  return(log.posterior)
}

### Setup
load(file = "CauchyData.RData")
data.cauchy = yVect

# a)
# prior
mu0 = 0
sigma0 = 10

# setup density sampling
width = 0.01
theta.grid = seq(0, 10, by=width)
postDensity = numeric()

# calculate logPosterior for each theta value
for(i in 1:length(theta.grid)) {
  postDensity[i] = logPostCauchyNormalPrior(theta.grid[i], gamma = 1, y = data.cauchy, mu0, sigma0)
}

# plot density
plot(theta.grid, (1/width)*exp(postDensity) / sum(exp(postDensity)), type="l",
     ylab="posterior density", xlab="theta")

# b)

logPostCauchy = function(theta.and.gamma, y, mu0, sigma.theta0, sigma.gamma0) {
  theta = theta.and.gamma[1]
  gamma = theta.and.gamma[2]
  
  # likelihood
  log.likelihood = sum(log(dCauchy(y, theta, gamma)))
  
  # prior theta
  log.prior.theta = dnorm(theta, mu0, sd=sigma.theta0, log = TRUE)
  
  # prior gamma
  log.prior.gamma = log(dlognormal(gamma, mu0, sigma.gamma0))
  
  # posterior
  log.posterior = log.likelihood + log.prior.theta + log.prior.gamma
  
  return(log.posterior)
}
theta.init = 0
gamma.init = 1
sigma.theta0 = 10
sigma.gamma0 = 1
optim.res = optim(c(theta.init, gamma.init), fn=logPostCauchy, 
                  data.cauchy, mu0, sigma.theta0, sigma.gamma0,
                  gr = NULL, method = "L-BFGS-B", lower=c(-Inf, 0.0001),
                  control=list(fnscale=-1), hessian=TRUE)

# c)
# extract mean and covariance matrix
post.mean = optim.res$par # Posterior mean vector
post.cov = -solve(optim.res$hessian) # Posterior covariance matrix

# sample
nSamples = 10000
post.samples = rmvnorm(nSamples, mean = post.mean, sigma = post.cov)
quant.99 = post.samples[,1] + post.samples[,2]*tan(pi*(0.99-1/2)) 
hist(quant.99, 100, freq=FALSE)
lines(density(quant.99))



