y = c(43,44,45,46,46.5,47.5)
theta = seq(0,100,1/1000)
posterior = rep(0,length(theta))
for ( i in 1 : length(theta)){
    posterior[i] = prod(1 / (1+(y-theta[i])^2 ))
}
posterior = posterior / sum(posterior)
plot(theta, posterior,type= "l", lwd =1, ylab="Density", main="posterior distribution")

# sample from posterior 
set.seed(33)
hist(sample(theta,size=100,replace = T, prob = posterior), main= "predictive histogram")

# predictive distribution:
par(mfcol=c(1,2))
set.seed(33)
obs = rep(0,1000)
sample = sample(theta, size = 1000, replace = T, prob = posterior)
for (i in 1:1000){
    obs[i] = rcauchy(1, location = sample[i])
}
hist(obs,breaks= 50,xlab="New observations", main="predictive y distribution")
## marginal distribution of y
## is it equvilant to conditional on mean of theta
set.seed(33)
sample2 = sample(theta,1000, replace = T, prob = posterior)
theta0 = mean(sample2)
obs2 = rcauchy(n=1000, location = theta0)
hist(obs2,breaks=50)

# different!!!!!