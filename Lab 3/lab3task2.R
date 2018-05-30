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
