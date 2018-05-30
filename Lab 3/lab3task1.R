#### Lab 3 Task 1: Normal model, mixture of normal model with semi-conjugate prior.
### Libraries

### Setup
# data
data.rain = read.table("rainfall.dat")$V1
n = dim(data.rain)[1]

# prior
mu.0 = mean(data.rain)
tau2.0 = 10
v.0 = 100 # << n so does not really matter
sigma2.0 = var(data.rain)

### Functions
rScaledInvChi2 = function(v, sigma2) {
  chi2 = rchisq(1, v)
  scaled.inv.chi2 = v * sigma2 / chi2
  return(scaled.inv.chi2)
}

rDirchlet = function(alpha, nAlloc) {
  n.categories = length(alpha)
  pi.draws = numeric()
  for(i in 1 : n.categories) {
    pi.draws[i] = rgamma(1, shape=alpha[i] + nAlloc[i], rate = 1) # shape = alpha, rate = beta
  }
  pi.draws = pi.draws / sum(pi.draws)
  return(pi.draws)
}

allocI = function(I) {
  n = dim(I)[1]
  alloc = numeric()
  for(i in 1:n) {
    alloc[i] = which(I[i,]==1)
  }
  return(alloc)
}

gibbsSamples = function(nSamples, x, mu.0, tau2.0, v.0, sigma2.0) {
  # set values
  mu.sample = mu.0
  sigma2.sample = sigma2.0
  n = length(x)
  
  # setup sampling
  mu.samples = numeric()
  sigma2.samples = numeric()
  
  # gibbs sampling
  v.n = v.0 + n
  for(i in 1:nSamples) {
    # sample sigma2
    sigma2.n = (v.0 * sigma2.0 + sum( (x - mu.sample)^2) ) / v.n
    sigma2.sample = rScaledInvChi2(v.n, sigma2.n)
    
    # sample mu
    w = (n / sigma2.sample) / (n / sigma2.sample + 1 / tau2.0)
    mu.n = w * mean(x) + (1 - w) * mu.0
    tau2.n = 1 / (n / sigma2.sample + 1 / tau2.0)
    mu.sample = rnorm(1, mean = mu.n, sd = sqrt(tau2.n))
    
    # save result
    sigma2.samples[i] = sigma2.sample
    mu.samples[i] = mu.sample
  }
  
  return(data.frame(mu = mu.samples, sigma2 = sigma2.samples))
}

### Implementation
# a) i)
nSamples = 4000
gibbs.samples = gibbsSamples(nSamples, data.rain, mu.0, tau2.0, v.0, sigma2.0)

# a) ii)
# posterior density
plot(gibbs.samples$mu, gibbs.samples$sigma2, col=rgb(0,0,0,0.3),
     xlab="mu", ylab="sigma2")

# trajectories
plot(gibbs.samples$mu, type="l", col="blue", xlab="sample", ylab="mu", main="trajectories for mu")
plot(gibbs.samples$sigma2, type="l", col="orange", xlab="sample", ylab="sigma2", main="trajectories for mu")

# estimate the paramateters of gibbs sampling using mean of posterior samples
mu.gibbs = mean(gibbs.samples$mu)
sigma2.gibbs = mean(gibbs.samples$sigma2)

# b)

# setup
x = as.matrix(data.rain)
nObs = length(x)
nComp = 2
nIterations = 100
prob.obs.comp = numeric()

# priors
alpha = 10 * rep(1, nComp)
muPrior = mu.0 * rep(1, nComp)
tau2Prior = tau2.0 * rep(1, nComp)
sigma2Prior = tau2.0 * rep(1, nComp)
vPrior = v.0 * rep(1, nComp)

# initial values
I = t(rmultinom(nObs, size = 1, prob = rep(1, nComp) / nComp))
muComp = numeric()
# muComp = quantile(x, probs = seq(0, 1, length = nComp))
# muComp = seq(min(x), max(x), length = nComp)
sigma2Comp = rep(var(x), nComp)

for(k in 1:nIterations) {
  alloc = allocI(I) # map allocation
  nAlloc = colSums(I) # calculate number of obs allocation to each component
  
  # update mu
  for(j in 1:nComp) {
    # w
    precData = nAlloc[j] / sigma2Comp[j]
    precPrior = 1 / tau2Prior[j]
    precPost = precData + precPrior
    w = precData / precPost 
    
    # mu.n
    muPost = w * mean(x[alloc==j]) + (1 - w) * muPrior[j] 
    
    # tau2.n
    tau2Post = 1 / precPost
    
    # mu
    muComp[j] = rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  # update sigma2
  for(j in 1:nComp) {
    # v.n
    vPost = vPrior[j] + nAlloc[j] 
    
    # sigma2.n
    sigma2Post = (vPrior[j] * sigma2Prior[j] + sum( (x[alloc==j] - muComp[j])^2) ) / vPost
    
    # sigma2
    sigma2Comp[j] = rScaledInvChi2(vPost, sigma2Post)
  }
  
  # update allocation
  pi = rDirchlet(alpha, nAlloc) # sample pi
  for(i in 1:nObs) {
    
    # calculate probability using sampled pi and density of given observation for corresponding component
    for(j in 1:nComp) {
      prob.obs.comp[j] = pi[j] * dnorm(x[i], mean = muComp[j], sd=sqrt(sigma2Comp[j]))
    }
    
    # randomize allocation for next iteration
    I[i,] = rmultinom(1, size = 1, prob = prob.obs.comp / sum(prob.obs.comp))
  }
}

# plot density curves
# plot(density(x), xlab="Precipation", main="Final fitted density") # plot density of data
hist(x, freq=FALSE, breaks=20, xlab="Precipation", main="Final fitted density")
x.grid = seq(-10, 400, 0.5)
densityComp1 = dnorm(x.grid, mean=muComp[1], sd=sqrt(sigma2Comp[1])) # density of comp1
densityComp2 = dnorm(x.grid, mean=muComp[2], sd=sqrt(sigma2Comp[2])) # density of comp2
lines(x.grid, pi[1]*densityComp1+pi[2]*densityComp2, col="red") # add curve of both densities, pi portioned
lines(x.grid, dnorm(x.grid, mean=mean(x), sd=sd(x)), col="blue") # add curve for a normal density estimation
legend("topright", legend=c("data histogram", "mixture density", "normal density"), col=c("black", "red", "blue"), lty=1)

# c) 
hist(x, freq=FALSE, breaks=20, ylim=c(0,0.035))
lines(x.grid, dnorm(x.grid, mean=mu.gibbs, sd=sqrt(sigma2.gibbs)), col="blue")
lines(x.grid, pi[1]*densityComp1 + pi[2]*densityComp2, col="red")
