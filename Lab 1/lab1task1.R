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
