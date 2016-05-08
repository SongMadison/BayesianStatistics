
bioassay <- data.frame(list(
  xi = c(-.86,-.30,-.05,.73),
  ni = c(5,5,5,5),
  yi = c(0,1,3,5)
  ))

model1 <- glm(cbind(yi,ni-yi)~xi, family =binomial(logit), data = bioassay)
summary(model1)

model2 <- glm(yi/ni~xi, weights = ni, family = binomial(logit), data = bioassay)
summary(model2)

prior <- function(alpha, beta) 
  exp(-2/3* (alpha^2/4+ (beta-10)^2/100- alpha*(beta-10)/20))

dens <- function(alpha,beta,y,n,x){
  p = 1
  for (i in 1:length(x)){
    t <- exp(alpha + beta*x[i])
    p = p* (t/(1+t))^y[i] * (1/(1+t))^(n[i]-y[i])
  }
  return (p)
}

post.unnorm <- function(alpha,beta,y=bioassay$yi,n= bioassay$ni,x=bioassay$xi) 
  prior(alpha,beta) * dens(alpha,beta,y,n,x)

alpha = seq (-5,10,length.out =1000); beta =seq(-10,40,length.out = 1000)
z = outer(alpha,beta, FUN=post.unnorm)
z <- z/(sum(z)*15/999*50/999)
contour(alpha,beta,z)

#simulate alpha, beta from posterior
alpha.mar <- rowSums(z)/sum(z)
alpha.draw <- sample(alpha,size = 1000,replace = T, prob = alpha.mar)
beta.draw <- numeric(1000)
for( i in 1:1000){
  beta.mar <- z[which(alpha == alpha.draw[i]),]/sum(z[which(alpha == alpha.draw[i]),])
  beta.draw[i] <- sample(beta,1,replace = T, prob = beta.mar)
}
alpha.draw <- jitter(alpha.draw,amount = 15/999)/2
beta.draw <- jitter(beta.draw, amount = 50/999/2) 
plot(alpha.draw,beta.draw,xlim = c(-5,10), ylim = c(-10,40) ,pch = ".")

sum(beta.draw>0)/1000

hist(-alpha.draw/beta.draw,breaks =50)

prior.m <- outer(alpha, beta, FUN = prior)
dat <- data.frame(list (alpha= rep(alpha,1000), beta= rep(beta, each= 1000),
                  post = as.vector(z), pri = as.vector(prior.m) ))
library(ggplot2)
ggplot(dat, aes(x= alpha, y =beta))+
  geom_contour(aes(z= post), col ="red", xlim= c(-5,10), ylim=c(-10,40))+
  geom_contour(aes(z= pri), col ="blue")
contour(alpha,beta,prior.m)





###########################################
###### Laplace approximation
#############################################################
#I = infomatrix(alpha, beta)
# N(MLE, I(MLE))

X <- matrix(c(rep(1,4),-0.86,-0.30,-.05,.73), ncol =2)
logit<-function(t) exp(t)/(1+exp(t))
beta.hat <- matrix(c(0.847,7.749))
W <- logit(X%*%beta.hat)
W <- W*(1-W)
I = t(X)%*%diag(as.vector(W))%*%X*5

solve(I)
### some mistake, not sure!!
dens1 <- function(alpha,beta,y,n,x){  
  t <- exp(alpha + beta*x)
  p = (t/(1+t))^y * (1/(1+t))^(n-y)
  return (p)
}
post.unnorm1 <- function(alpha,beta,y=bioassay$yi,n= bioassay$ni,x=bioassay$xi) 
  prior(alpha,beta) * prod(dens1(alpha,beta,y,n,x))
z <- function(alpha,beta){
  mat <- matrix(nrow = length(alpha), ncol =length(beta))
  for ( i in 1:length(alpha))
    for ( j in 1 :length(beta))
      mat[i,j]  <- post.unnorm(alpha[i],beta[j])
  return(mat)
}
alpha = seq (-2,2,0.005); beta =seq(-4,4,0.005)
mat <- z(alpha,beta)
contour(alpha,beta, z=mat)
