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

