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
