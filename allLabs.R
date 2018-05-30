α = alpha
β = beta
Γ = gamma
ε = epsilon
θ = theta
λ = lambda
μ = mu
Φ = phi
ψ = psi
Ω = Omega
ω = omega
χ = chi
Σ = Sigma
σ = sigma
########################## Lab1 ##########################
### Meta-info
## Beta distribution
## Simulations
## Gini coefficient
## Credible interval
## Highest Posterior Density HPD

# Task 1: Bernoulli ... again
# a) Draw random numbers from Beta distribution and graphical verification
# of posterior
# b) Simulate to compute posterior prob of Pr(theta < 0.4)
# c) Compute log-oods posterior distribution

# Task 2: Log-normal distribution and the Gini coefficient
# a) Simulate 1000 draws from posterior or theta2. Compare with real value
# b) Compute posterior distribution of Gini coefficient G
# c) Compute 95% equal tail credible interval of Gini coefficient G.
# Doing a kernal density estimate
# Compute 95% Highest Posterior Density interval (HPD) of G

# Task 3: Bayesian inference for the concentration parameter in the von Mises distributio
# a) Plot posterior distribution of kappa for wind direction data
# b) Find approximate posterior mode of kappa

#####################################################################################################

#### Lab 1 Task 1: Bernoulli... Again

### Setup
s = 14
n = 20
f = n - s # = 6
alpha0 = 2
beta0 = 2

### Functions
logOdds = function(theta) {
  return(log(theta/(1-theta)))
}

### Implementation
## a) Draw samples from posterior theta and verify graphically the convergence as sample size grows

# Setup for sampling
nDraws = 100000
nSeq = seq(10, nDraws, 100)
thetaDraw.means = numeric()
thetaDraw.sds = numeric()

# sample thetas, store mean and standard deviations
for(i in 1:length(nSeq)) {
  # sample n thetas
  thetaDraw = rbeta(nSeq[i], alpha0 + s, beta0 + f)
  
  # save mean of samples
  thetaDraw.means[i] = mean(thetaDraw)
  
  # save standard deviation of samples
  thetaDraw.sds[i] = sd(thetaDraw)
}

# plot posterior mean and standard deviation for rising n's to witness convergence
# mean
true.mean = (s + alpha0) / (s + alpha0 + f + beta0)
plot(nSeq, thetaDraw.means, type="l", col="blue")
abline(h=true.mean, col="green")

# standard deviation
post.alpha = s + alpha0
post.beta = f + beta0
true.sd = sqrt((post.alpha * post.beta) / ((post.alpha + post.beta)^2 * (post.alpha + post.beta + 1)))
plot(nSeq, thetaDraw.sds, type="l", col="red")
abline(h=true.sd, col="black")

## b) Compute posterior probability P(theta < 0.4) using nDraws = 10000
nDraws = 10000
theta.samples = rbeta(nDraws, alpha0 + s, beta0 + f)
sampled.prob = sum(ifelse(theta.samples < 0.4, 1, 0)) / nDraws # 0.0033
actual.prob = pbeta(0.4, alpha0 + s, beta0 + f) # 0.003972681

## c) Compute posterior density of the log-odds = phi = log(theta / 1 - theta) using sampling
nDraws = 10000
theta.samples = rbeta(nDraws, alpha0 + s, beta0 + f)
phiSamples = logOdds(theta.samples)
plot(density(phiSamples))

#### Lab 1 Task 2: Log-normal distribution and the Gini coefficient

### Setup
y = c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
n = length(y)
mu = 3.5

### Functions

### Implementation
## a) draw 10000 samples from Scaled-Inv Chi2(tau2, n)
tau2 = sum((log(y) - mu)^2) / n
nDraws = 10000
chi2.samples = rchisq(nDraws, n)
scaled.chi2.samples = n * tau2 / chi2.samples
sigma2.samples = scaled.chi2.samples

# plot posterior density
plot(density(sigma2.samples))

# b) sample Gini Coefficients
sigma.samples = sqrt(sigma2.samples)
G = 2 * pnorm(sigma.samples / sqrt(2)) - 1

# c) Credible interval and HPD for Gini
# calculate credible interval
G.sorted = sort(G)
CI.lower = G.sorted[length(G.sorted)*0.025]
CI.upper = G.sorted[length(G.sorted)*0.975]

# plot credible interval boundries
plot(density(G))
abline(v=CI.lower, col="green")
abline(v=CI.upper, col="green")

# calculate highest posterior density interval
G.density = density(G)
G.density = data.frame(x=G.density$x, y=G.density$y)
G.density.sorted = G.density[order(G.density$y),]
G.percentages = cumsum(G.density.sorted$y) / sum(G.density.sorted$y)
G.HPD = G.density.sorted[which(G.percentages >= 0.05),] # extract the top 95% density
HPD.lower = min(G.HPD$x)
HPD.upper = max(G.HPD$x)

# plot HPD interval boundries
plot(density(G))
abline(v=HPD.lower, col="blue")
abline(v=HPD.upper, col="blue")

#### Lab 1 Task 3: Bayesian inference for the concentration parameter in the von Mises distribution

### Setup
y = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
n = length(y)
mu = 2.39

### Functions
likelihoodVonMises = function(k, y, mu) {
  return(prod(exp(k * cos(y - mu)) / (2 * pi * besselI(k, 0))))
}

### Implementation
# setup grid
k.grid = seq(0, 10, 0.01)

# calculate likelihoods
likelihoods = numeric() 
for(i in 1:length(k.grid)) {
  likelihoods[i] = likelihoodVonMises(k.grid[i], y, mu)
}
# calculate prior of kappa with lambda = 1
# prior = dexp(k.grid) # rate = 1 is default
prior = exp(-k.grid)

# calculate posterior
posterior = likelihoods * prior

# plot posterior
plot(k.grid, posterior, type="l")

# k corresponding to max posterior density
k.max = k.grid[which.max(posterior)]

########################## Lab 2 ##########################
## Meta-info
# Linear regression
# Polynomial regression
# Logistic regressions
# Credible interval
# Maximum likelihood
# Optim
# Hessian
# Mode of beta
# Predictive distributions

# Task 1: Linear and polynomial regression
# a) Set the prior hyperparameters µ0, Ω0, ν0 and σ2 to sensible values
# b) Check if your prior from a) is sensible
# Simulate draws from joint prior and compute regression curve of each draw
# c) Simulates from the joint posterior distribution of β0, β1,β2 and σ2
# Plot: 
# • Posterior mean of simulations
# • Curve of lower and upper credible interval of f(time)
# d) Simulate highest expected temperatures
# Simulate from posterior distribution of time with highest expected temperatures
# What to do to mitigate risk of overfitting high order polynomial regression?

# Task 2: Posterior approximation for classification with logistic regression
# a) Fit logistic regression using maximum likelihood estimations
# b) Approximate the posterior distribution of the 8-dim parameter vector β with a multivariate normal distribution
# c) Simulates from the predictive distribution of the response variable in a logistic regression

#####################################################################################################

#### Lab 2 Task 1: Linear and polynomial regression
### Libraries
library(mvtnorm)

### Setup
data.temp = read.table("TempLinkoping.txt", header=TRUE)
time = data.temp$time
X = cbind(1, data.temp$time, data.temp$time^2)
y = data.temp$temp
n = dim(data.temp)[1]
p = dim(X)[2]
# a)
beta.0 = c(-4, 110, -110)
Omega.0 = diag(c(2, 2, 2))
v.0 = 100
sigma2.0 = 9

### Functions
rscaledinvchisq = function(n, v, sigma2) {
  
  # draw from chi2(v)
  chi2 = rchisq(n, v) 
  
  # draw from scaled inv-chi2(v, sigma2)
  scaled.inv.chi2 = v * sigma2 / chi2
  
  return(scaled.inv.chi2)
}

### Implementation
# b)
# new parameters to tweak
mu.0 = c(-4, 100, -100)
# sigma0 changes the shape drastically:
# large values -> small spread
# small values -> large spread
Omega.0 = diag(c(0.5, 0.5, 0.5))
v.0 = 100
sigma2.0 = 4

# start empty plot
plot(NULL, xlim=c(0,1), ylim=c(-10,30), ylab="temperature", xlab="time")

# number of samples
nSamples = 1000

# sample sigma2s
sigma2.samples = rscaledinvchisq(nSamples, v.0, sigma2.0)

for(sigma2.sample in sigma2.samples) {
  # sample beta
  beta.sample = t(rmvnorm(1, mean = mu.0, sigma = sigma2.sample * solve(Omega.0)))
  
  # add regression curve to plot
  lines(time, X%*%beta.sample, col=rgb(0,0,0,0.2))
}

# c)
# calculate posterior parameters
beta.hat = solve(t(X) %*% X) %*% t(X) %*% y
mu.n = solve(t(X) %*% X + Omega.0) %*% (t(X) %*% X %*% beta.hat + Omega.0 %*% mu.0)
Omega.n = t(X) %*% X + Omega.0
v.n = v.0 + n
sigma2.n = ((v.0 * sigma2.0 + t(y) %*% y + t(mu.0) %*% Omega.0 %*% mu.0 - t(mu.n) %*% Omega.n %*% mu.n) / v.n)[1]

# setup sampling
plot(data.temp)
beta.samples = matrix(NA, nSamples, p)
y.samples = matrix(NA, nSamples, n)

# sample sigma2s
sigma2.samples = rscaledinvchisq(nSamples, v.n, sigma2.n)

for(i in 1:nSamples) {
  # sample beta
  beta.samples[i,] = rmvnorm(1, mean = mu.n, sigma = sigma2.samples[i] * solve(Omega.n))
  
  # calculate y for sample
  y.samples[i,] = X %*% beta.samples[i,]
  
  # add posterior regression to data plot
  lines(time, y.samples[i,], col=rgb(0,0,0,0.2))
}

# extract credible interval
interval = function(x) {
  quantile(x, probs=c(0.05, 0.95), na.rm=T)
}
post.CI = apply(y.samples, 2, interval)
post.mean = apply(y.samples, 2, mean)

# plot data with credible interval
plot(data.temp)
lines(time, post.CI[1,], col="red")
lines(time, post.CI[2,], col="red")
lines(time, post.mean, col="green")

# d)
# setup sampling
y.best.samples = numeric();
time.best.samples = numeric()

# sample from time with highest temperature
for(i in 1:nSamples) {
  # extract sample
  beta.sample = beta.samples[i,]
  
  # calculate best time
  time.best.samples[i] = -beta.sample[2]/(2*beta.sample[3])
  
  # calculate and save temperature for best time
  y.best.samples[i] = c(1, time.best.samples[i], time.best.samples[i]^2) %*% beta.sample
}

# hist temperature of best time
hist(y.best.samples, breaks=20) # around 15.8
hist(time.best.samples*366, breaks=20) # day 183

# e)
# To mitigate the problem of overfitting when using high order polynomial models, 
# it’s reasonable to replace the prior with a ridge regression, 
# where mu=0 and Omega0 is a diagonal matrix with lambda values. 
# To limit the higher orders, larger values for lambdas 
# corresponding to the higher orders should be chosen:

#### Lab 2 Task 2: Posterior approximation for classification with logistic regression
### Libraries
library(mvtnorm)

### Setup
# data
data.work = read.table("WomenWork.dat", header=TRUE)
X = as.matrix(data.work[,-1])
y = data.work[,1]
n = length(y)
p = dim(X)[2]

# priors
mu.0 = rep(0, p)
tau = 10
Sigma.0 = tau^2 * diag(p)

### Functions
logPostLogistic = function(beta, X, y, mu.prior, Sigma.prior) {
  
  # log likelihood
  log.likelihood = sum( y * X %*% beta - log(1 + exp(X %*% beta)) )
  if (abs(log.likelihood) == Inf) log.likelihood = -20000
  
  # log prior
  log.prior = dmvnorm(beta, mean = mu.prior, sigma = Sigma.prior, log = TRUE)
  
  # log posterior
  log.posterior = log.likelihood + log.prior
  
  return(log.posterior)
}

### Implementation
# a) 
logistic.model = glm(Work ~ 0 + ., data = data.work, family = binomial)
y.pred = ifelse(logistic.model$fitted.values > 0.5, 1, 0)
mcr = sum(abs(y.pred-data.work$Work)) / n

# b) 
optim.res = optim(mu.0, # inital values of beta
                  logPostLogistic, # function to be minimized (maximized)
                  X, y, mu.0, Sigma.0, # input parameters for the function
                  gr=NULL, # function to return the gradiant to
                  method="BFGS", # some gradiant function
                  control=list(fnscale=-1), # scale function result with -1 to maximize
                  hessian=TRUE) # if return hessian

# extract results
beta.mode = optim.res$par
beta.inv.neg.hessian = -solve(optim.res$hessian)

## CI by sampling
# sample betas
nSamples = 1000
beta.samples = rmvnorm(nSamples, mean = beta.mode, sigma = beta.inv.neg.hessian)

# calculate credible.intervals
credible.intervals = as.data.frame(apply(beta.samples, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
colnames(credible.intervals) = names(data.work[,-1])
print(credible.intervals)
#        Constant  HusbandInc EducYears   ExpYears  ExpYears2         Age NSmallChild  NBigChild
# 2.5%  -2.353987 -0.04783091 0.0287374 0.03394872 -0.6113804 -0.13420979  -2.1462354 -0.2941665
# 97.5%  3.635838  0.01095920 0.3295878 0.29935404  0.3237674 -0.02852173  -0.5726629  0.2402660
# 0 is not in the CI for NSmallChild -> it's an important determinant of the prob a woman works

## CI by calculation
beta.NSmallChild.mean = beta.mode[7]
beta.NSmallChild.sd = sqrt(beta.inv.neg.hessian[7,7])
CI.lower.NSmallChild = beta.NSmallChild.mean - 1.96 * beta.NSmallChild.sd # -2.121425
CI.upper.NSmallChild = beta.NSmallChild.mean + 1.96 * beta.NSmallChild.sd # -0.5968414

# c)
woman = c(1, 10, 8, 10, (10/10)^2, 40, 1, 1)

# sigmoid function
sigmoid = function(x) {
  return(exp(x) / (1 + exp(x)))
}

# returns predictions of logisitic regression
glmPred = function(n, x, mu, Sigma, linkFun) {
  beta.samples = rmvnorm(n, mean = mu, sigma = Sigma)
  lin.reg.samples = x %*% t(beta.samples)
  prob.pred = linkFun(lin.reg.samples)
  return(prob.pred)
}

work.probs = glmPred(nSamples, woman, beta.mode, beta.inv.neg.hessian, sigmoid)
work.pred = ifelse(work.probs > 0.5, 1, 0)
hist(work.pred)

########################## Lab 3 ##########################
## Meta-info
# Gibbs sampling
# MCMC
# Mixure normal models
# STAN
# Credible interval
# Effective samples

# Task 1: Normal model, mixture of normal model with semi-conjugate prior
# a) Implement Gibbs sampler that simulate from a joint posterior. Data is normal
# Evaluate Gibbs model graphically
# b) Implement Gibbs sampler that simulates form a joint posterior when mixture normal models
# Evaluate convergance of model

# Task 2: Time series models in Stan
# a) Write a function in R that simulates data from given AR(1)-process
# b) Treat µ, φ and σ2 as unknowns and estimate them using MCMC. 
# Report credible interval and effective samples
# c) Same as a) and b) but with new data
# d) Change σ2 to formative prior

##########################################################

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

#### Lab 3 Task 2: Time series models in Stan
### Libraries
library(rstan)

### Setup
mu = 10
sigma2 = 2
T = 200
phis = seq(-1, 1, length = 41)
nSamples = length(phis)

### Functions
ARSimulation = function(phi, mu, sigma2, T) {
  x = numeric()
  x[1] = mu
  for(t in 2:T){
    epsilon = rnorm(1, mean = 0, sd = sqrt(sigma2))
    x[t] = mu + phi * (x[t-1] - mu) + epsilon
  }
  return(x)
}

### Implementation
# a)
AR.samples = matrix(NA, nSamples, T)
for(i in 1:nSamples) {
  AR.samples[i,] = ARSimulation(phis[i], mu, sigma2, T)
}

# plot results
par(mfrow=c(2,2))
plot(AR.samples[1,], xlab="t", ylab="x_t", main=paste("phi =", phis[1]), type="l")
plot(AR.samples[15,], xlab="t", ylab="x_t", main=paste("phi =", phis[15]), type="l")
plot(AR.samples[27,], xlab="t", ylab="x_t", main=paste("phi =", phis[27]), type="l")
plot(AR.samples[40,], xlab="t", ylab="x_t", main=paste("phi =", phis[40]), type="l")
par(mfrow=c(1,1))

# b)
# simulate samples
AR.0.3 = ARSimulation(0.3, mu, sigma2, T)
AR.0.95 = ARSimulation(0.95, mu, sigma2, T)

# stan model
ARStanModel = 
  'data  {
    int<lower=0> T;
    vector[T] x;
}
parameters {
    real mu;
    real sigma2;
    real phi;
}
model {
    for(t in 2:T) {
      x[t] ~ normal(mu + phi * (x[t-1] - mu), sqrt(sigma2));
    }
}'

# fit stan model
fit.0.3 = stan(model_code = ARStanModel, data = list(x = AR.0.3, T = T))
fit.0.95 = stan(model_code = ARStanModel, data = list(x = AR.0.95, T = T))

# save summaries
summary.0.3 = summary(fit.0.3)$summary
summary.0.95 = summary(fit.0.95)$summary

# plot interesting info about stan fit
# traceplot(fit) # traceplot of all chains
# pairs(fit) # bivariate posterior plots
# save plot
# png('myPlot.png')
# plot(rnorm(4)) # Change this to the plot commands you want 
# dev.off()

# i) extract mean and confidence intervals (mean lower_bound upper_bound n_eff)
mu.info.0.3 = c(summary.0.3["mu", "mean"], summary.0.3["mu", "2.5%"], 
                summary.0.3["mu", "97.5%"], summary.0.3["mu", "n_eff"])
# 9.985365  9.705207 10.268224 (compared to 10) 3666.931649

sigma2.info.0.3 = c(summary.0.3["sigma2", "mean"], summary.0.3["sigma2", "2.5%"], 
                    summary.0.3["sigma2", "97.5%"], summary.0.3["sigma2", "n_eff"])
# 2.038662 1.665407 2.504641 (compared to 2) 2947.096362

phi.info.0.3 = c(summary.0.3["phi", "mean"], summary.0.3["phi", "2.5%"], 
                 summary.0.3["phi", "97.5%"], summary.0.3["phi", "n_eff"])
# 0.2585258 0.1195351 0.3975016 (compared to 0.3) 4000.0000000

mu.info.0.95 = c(summary.0.95["mu", "mean"], summary.0.95["mu", "2.5%"], 
                 summary.0.95["mu", "97.5%"], summary.0.95["mu", "n_eff"])
# 8.785386  1.461110 13.539001 (compared to 10) 458.824080

sigma2.info.0.95 = c(summary.0.95["sigma2", "mean"], summary.0.95["sigma2", "2.5%"], 
                     summary.0.95["sigma2", "97.5%"], summary.0.95["sigma2", "n_eff"])
# 1.816415 1.490158 2.232466 (compared to 2) 2311.436439

phi.info.0.95 = c(summary.0.95["phi", "mean"], summary.0.95["phi", "2.5%"], 
                  summary.0.95["phi", "97.5%"], summary.0.95["phi", "n_eff"])
# 0.9301223 0.8664641 0.9952547 (compared to 0.95) 1035.3806965

# ii) joint posterior of phi and mu
samples.0.3 = extract(fit.0.3)
samples.0.95 = extract(fit.0.95)

# extract mu and phi
mu.samples.0.3 = samples.0.3$mu
phi.samples.0.3 = samples.0.3$phi
mu.samples.0.95 = samples.0.95$mu
phi.samples.0.95 = samples.0.95$phi

# plot joint posteriors
plot(mu.samples.0.3, phi.samples.0.3, col=rgb(0,0,0,0.25))
plot(mu.samples.0.95, phi.samples.0.95, col=rgb(0,0,0,0.25))

# c) 
# setup
data.campy = read.table("campy.dat", header=TRUE)[,1]
n = length(data.campy)

PoissonARStanModel = 
  '
data {
    int<lower=0> N;
    int c[N];
}
parameters {
    real mu;
    real sigma2;
    real phi;
    vector[N] x;
}
model {
    for(t in 2:N) {
      x[t] ~ normal(mu + phi * (x[t-1] - mu), sqrt(sigma2));
      c[t] ~ poisson(exp(x[t]));
    }
}
'
# fit model
poissonAR.fit = stan(model_code=PoissonARStanModel, data=list(N=n, c=data.campy))

# save summary
poissonAR.summary = summary(poissonAR.fit)$summary

# extract upper bound, lower bound and posterior mean
x.upper = poissonAR.summary[4:n, "97.5%"]
x.lower = poissonAR.summary[4:n, "2.5%"]
x.mean = poissonAR.summary[4:n, "mean"]

# plot data with confidence intervals
plot(data.campy, xlab="t", ylab="number of cases of campylobacter")
lines(exp(x.upper), col="red")
lines(exp(x.lower), col="red")
lines(exp(x.mean), col="green")

# d)
SmoothPoissonARStanModel = 
  '
data {
    int<lower=0> N;
    int c[N];
}
parameters {
    real mu;
    real sigma2;
    real phi;
    vector[N] x;
}
model {
    phi ~ uniform(-1, 1);
    sigma2 ~ scaled_inv_chi_square(1000, 0.02);
    for(t in 2:N) {
        x[t] ~ normal(mu + phi * (x[t-1] - mu), sqrt(sigma2));
        c[t] ~ poisson(exp(x[t]));
    }
}
'
# settings for rstan
nBurnin = 1000
nSamples = 1000
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# fit model
smoothPoissonARFit = stan(
  model_code=SmoothPoissonARStanModel, 
  data=list(N=n, c=data.campy),
  warmup = nBurnin,
  iter=(nBurnin+nSamples),
  chains=2
)

# extract summary
smooth.summary = summary(smoothPoissonARFit)$summary

# extract upper bound, lower bound and posterior mean
smooth.x.upper = smooth.summary[4:n, "97.5%"]
smooth.x.lower = smooth.summary[4:n, "2.5%"]
smooth.x.mean = smooth.summary[4:n, "mean"]

# plot data with confidence intervals
plot(data.campy, xlab="t", ylab="number of cases of campylobacter")
lines(exp(smooth.x.upper), col="red")
lines(exp(smooth.x.lower), col="red")
lines(exp(smooth.x.mean), col="green")

########################## Lab 4 ##########################
## MCMC
## Metropolis algorithm
## Poisson regression
## Hessian
## Mode
## Maximum likelihood estimator
## Program general functions

# Task 1: Poisson regression - the MCMC way
# a) Obtain the maximum likelihood estimator of β in the Poisson regression model
# Find significant covariates
# b) Bayesian analysis of the Poisson regression
# Find mode and hessian of Beta
# c) Simulate from the actual posterior of β using the Metropolis algorithm and 
# compare with the approximate results in b)
# Program general function
# d) Use MCMC draws from c) to simulate from the predictive distribution of
# the number of bidders in a new auction

#####################################################################################################

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
