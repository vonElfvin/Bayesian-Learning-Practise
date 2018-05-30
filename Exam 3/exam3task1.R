#### Exam 3 Task 1: Bayesian inference for proportions data
### Setup
load(file = 'yProportions.RData')
data.prop = yProp
n = length(data.prop)

# a)
logPostBetaExpPrior = function(theta, y) {
  # likelihood
  log.likelihood = sum(dbeta(y, theta, theta, log = TRUE))
  
  # prior
  log.prior = dexp(theta, rate = 1, log = TRUE)
  
  # posterior
  log.posterior = log.likelihood + log.prior
  
  return(log.posterior)
}

theta.grid = seq(0.01, 15, length = 1000)
logPost = numeric()
for(i in 1:length(theta.grid)) {
  logPost[i] = logPostBetaExpPrior(theta.grid[i], data.prop)
}

# calculate width to scale density to integrate to 1 over theta
width = theta.grid[2] - theta.grid[1]

# plot density
plot(theta.grid, exp(logPost) / (sum(exp(logPost)) * width), type = "l")

# zero-one loss -> posterior mode is the optimal Bayes point esitmator
post.mode = theta.grid[which.max(logPost)] # 4.481491

# b)
logPostBetaExpPriors = function(params, y) {
  # extract thetas
  theta1 = params[1]
  theta2 = params[2]
  
  # likelihood
  log.likelihood = sum(dbeta(y, theta1, theta2, log = TRUE))
  
  # prior
  log.prior.theta1 = dexp(theta1, rate = 1, log = TRUE)
  log.prior.theta2 = dexp(theta2, rate = 1, log = TRUE)
  
  # posterior
  log.posterior = log.likelihood + log.prior.theta1 + log.prior.theta2
  
  return(log.posterior)
}

optim.res = optim(c(1,1), logPostBetaExpPriors, data.prop,
                  gr = NULL, method = "L-BFGS-B", lower = c(0.001, 0.001),
                  control = list(fnscale = -1), hessian = TRUE)
post.mean = optim.res$par
post.covariance = -solve(optim.res$hessian)
