### Exam 3 Task 4: Prediction and decision
### Setup
y = c(195, 191, 196, 197, 189)
n = length(y)
sigma = 10
nSamples = 1000

# a)

# sample max weight for next day
# mu pred = mean(data)
# sd pred = sigma^2 * (1 + 1/n)
maxWeightNextDay.samples = numeric()
for(i in 1:nSamples) {
  maxWeightNextDay.samples[i] = rnorm(1, mean = mean(y), sd = sqrt(sigma^2 * (1 + 1/n)))
}
hist(maxWeightNextDay.samples, 20)

# b)

maxWeightNextYear.samples = numeric()
for(i in 1:nSamples) {
  maxWeightNextYear.samples[i] = max(rnorm(365, mean = mean(y), sd = sqrt(sigma^2 * (1 + 1/n))))
}
hist(maxWeightNextYear.samples, 20)
prob.230 = sum(maxWeightNextYear.samples > 230) / nSamples # 0.149

# c)

lossFunc = function(a, y.max) {
  # the bridge does not fall: loss = cost of building
  if(y.max <= 10 * a) {
    return(a)
    
  # the bridge falls: loss = cost of building + 100
  } else {
    return(a + 100)
  }
}

a.grid = seq(10, 30, by=0.01)
nActions = length(a.grid)

action.loss.samples = matrix(NA, nActions, nSamples)
for(i in 1:nActions) {
  for(j in 1:nSamples) {
    action.loss.samples[i,j] = lossFunc(a.grid[i], maxWeightNextYear.samples[j])
  }
}
action.loss.post.mean = apply(action.loss.samples, 1, mean)
a.grid[which.min(action.loss.post.mean)]
plot(a.grid, action.loss.post.mean, type="l")
