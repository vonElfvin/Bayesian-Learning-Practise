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
     