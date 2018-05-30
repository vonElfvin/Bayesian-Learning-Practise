#### Exam 1 Task 1: Bayesian inference for the rice distribution

### Help Functions

# Random number generator for the Rice distribution
rrice = function(n = 1, theta = 1, psi = 1){
  x <- rnorm(n = n, mean = 0, sd = sqrt(psi))
  y <- rnorm(n = n, mean = theta, sd = sqrt(psi))
  return(sqrt(x^2+y^2))
}

### Setup
psi = 1
data.rice = riceData
x = data.rice

### Implementation
# a)

# Density function for rice distribution
logdrice = function(theta, x, psi) {
  logd = log(x / psi) - (x^2 + theta^2) / (2 * psi) + log(besselI(x * theta / psi, 0))
  return(logd)
}

# Log post density for the rice posterior
logPostRice = function(theta, psi, x) {
  # likelihood
  log.likelihood = sum(logdrice(theta, x, psi))
  
  # prior
  log.prior = 0 # prior is assumed to be constant
  
  # posterior
  log.posterior = log.likelihood + log.prior
  return(log.posterior)
}

# setup grid
grid.width = 0.01
theta.grid = seq(0.01, 3, by = grid.width)
logPostRice.grid = numeric()

# calculate log of the posterior of the rice distribution
for(i in 1:length(theta.grid)) {
  logPostRice.grid[i] = logPostRice(theta.grid[i], psi, x)
}
# plot the actual density (note grid.width and divided by sum)
plot(theta.grid, (1 / grid.width) * exp(logPostRice.grid) / sum(exp(logPostRice.grid)), type="l")

# b)
# init value (any would suffice)
theta.init = mean(x)

# optimize over log of the posterior
optim.res = optim(theta.init, logPostRice, psi, x,
                  gr = NULL, method = "L-BFGS-B", lower = 0, control = list(fnscale=-1), hessian=T)

# extract results
theta.mode = optim.res$par
theta.sd = sqrt(-1/optim.res$hessian[1,1])

# compare actual density to normal density approiximation
plot(theta.grid, (1 / grid.width) * exp(logPostRice.grid) / sum(exp(logPostRice.grid)), type="l")
lines(theta.grid, dnorm(theta.grid, mean=theta.mode, sd=sqrt(theta.sd)))
# Comment: decent result, longer tails than the original but could see practical use

# c)
# number of samples
nSamples = 10000

# sample p(theta|x)
theta.samples = rnorm(nSamples, mean = theta.mode, sd = sqrt(theta.sd))

# sample p(x[n+1] | theta[i])
x.samples = numeric()
for(i in 1:nSamples) {
  x.samples[i] = rrice(n=1, theta=theta.samples[i], psi=1)
}
# plot the distribution of the sampled x[t+1]
hist(x.samples, 40, freq=FALSE)
