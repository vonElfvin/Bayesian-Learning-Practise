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
  

