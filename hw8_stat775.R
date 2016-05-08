#ex2 : # bayesian one-way anova. 

## Modeling 
#  y_ij   = theta_i  + \epsion_ij
#theta_i ~ N(mu, sigma_t2),   \epsilon_ij ~ N(0, sigma.e2)
#prior:
# mu \sim N( mu, sigma.02)
# sigma.t2 ~ IG(a_1, b_1) , sigma.e2 ~ IG(a_2, b_2) ,,, IG is inverse-Gamma distribution
# hyperprior:
# mu ~ N(mu0, sigma.02)


rm( list = ls() )
#install.packages("MCMCpack")
library(MCMCpack)

y =c(7.298, 3.846,  2.434, 9.566, 7.990,
     5.220, 6.556,  0.608,11.788, -.892,
     0.110,10.386, 13.434, 5.510, 8.166,
     2.212, 4.852,  7.092, 9.288, 4.980,
     0.282, 9.014,  4.458, 9.446, 7.198,
     1.722, 4.782,  8.106, 0.758, 3.758
)
y <- matrix(y,nrow=6, byrow = T)
y.bar <- apply(y, MARGIN = 1, mean)



# prior parameters -- prior distribution
a_1 = 1; b_1 = 1   # sigma.t2
a_2 = 0; b_2 = 0   # sigma.e2
mu0 <- 0 ; sigma.02 <-  10^12
K=dim(y)[1]; J = dim(y)[2]
iterN = 75000 


## initialzing ...
mu <- 0.1
sigma.e2 = 1
sigma.t2 = 1
theta <- rnorm(K, mu, sigma.t2^0.5)
result <- matrix(0,nrow = 3+K, ncol = iterN)

result[,1] <- c(mu, sigma.e2, sigma.t2, theta)

#Gibbs sampling 

set.seed(775)
for( j in 2:iterN){
    ## mu
    mean.mu <- (sigma.t2*mu0+sigma.02*sum(theta))/(sigma.t2+K*sigma.02)
    sigma2.mu <- (sigma.t2*sigma.02) / (sigma.t2+K*sigma.02)
    mu <- rnorm(1, mean.mu, sigma2.mu^(0.5))
    
    ## theta_i
    mean.theta <- numeric(K); sigma2.theta <- numeric(K)
    for ( i in 1:K){
        mean.theta[i] <- J*sigma.t2/(J*sigma.t2+sigma.e2)* y.bar[i] 
                               + sigma.e2/(J*sigma.t2+sigma.e2)* mu
        sigma2.theta[i]<- sigma.t2*sigma.e2/(J*sigma.t2+sigma.e2)
        theta[i] <- rnorm(1, mean.theta[i], sigma2.theta[i]^(1/2))
    }
   
    ## sigma.e2
    ss <- 0
    for(i in 1:K){
        ss <- ss + sum((y[i,]-theta[i])^2)
    }
    ss <- ss/2
    sigma.e2 <- rinvgamma(n=1, shape = a_2 + K*J/2, scale = b_2 +ss) 
    
    ## sigma.t2
    sigma.t2 <- rinvgamma(n=1, shape =  a_1 + K/2, scale = b_1 + sum((theta-mu)^2)/2)
    
    result[,j] <- c(mu, sigma.e2, sigma.t2, theta)
}


result <- t(result)
result <- data.frame(result)
names(result) <- c("mu","sigma.e2","sigma.t2","theta.1","theta.2","theta.3","theta.4","theta.5",
                   "theta.6")
apply(  result[70000:75000,],MARGIN = 2, 
      function(x) quantile(x, probs = c(0.025,0.25,0.5,0.75,0.975)) )

plot(70000:75000,result[["theta.1"]][70000:75000], ylab = "theta", pch =".", type = "b")
