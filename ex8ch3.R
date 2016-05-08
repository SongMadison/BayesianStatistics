#exercise 8 chapter 3

y = c(16/58,9/90,10/48,13/57,19/103,20/57,18/86,17/112,35/273,55/64)
z = c(12/113,1/18,2/14,4/44,9/208,7/67,9/29,8/154)
y = y/(1+y)
z = z/(1+z)

theta.y =function(theta){
  p = 1
  for ( i in 1:10){
    p = p*gamma(10+theta)/gamma(theta)*y[i]^9*(1-y[i])^(theta-1)
  }
  return(p)
}

theta.z =function(theta){
  p = 1
  for ( i in 1:8){
    p = p*gamma(10+theta)/gamma(theta)*z[i]^9*(1-z[i])^(theta-1)
  }
  return(p)
}

ty = 1:100
tz = 1:120
c1 = sum(theta.y(ty))
py = theta.y(ty)/c1

c2 = sum(theta.z(tz))
pz = theta.z(tz)/c2

contour(ty,tz,outer(py,pz),xlab= "theta.y", ylab= "theta.z")
str(outer(py,pz))
?contour

## sampling from posterior distribution 
sample.y <- sample(1:100,1000, replace = T, prob = py)
sample.z <- sample (1:120, 1000, replace = T, prob  = pz)
diff.mu = 10/(sample.y+10) - 10/(sample.z+10)
hist(diff.mu, breaks =100)

####################################################################################################
# How to simulate from join distribution
##################################################################################3
if (!require("tikzDevice")){   
  install.packages("tikzDevice") 
  require("tikzDevice")
}
library(tikzDevice)

prior <- function(theta) {  # theta is a vector
  exp(-0.01 * sum(theta)) * 1e-4
}
dens <- function(x, theta) 
  (beta(theta[1],theta[2])^(-length(x)))*
      (prod(x)^(theta[1]-1)* prod(1-x)^(theta[2]-1))
post.unnorm <- function(x,theta)  prior(theta)*dens(x,theta)

m <- function(x, a, b){
  mat <- matrix (nrow = length(a), ncol= length(b))
  for ( i in 1:length(a))
    for (j in 1:length(b))
      mat[i,j] <- post.unnorm(x,c(a[i],b[j]))
  return(mat)
}

y = c(16/58,9/90,10/48,13/57,19/103,20/57,18/86,17/112,35/273,55/64)
z = c(12/113,1/18,2/14,4/44,9/208,7/67,9/29,8/154)
y = y/(1+y)
z = z/(1+z)

theta.y.1 <- seq (1,8,0.005); theta.y.2 <- seq(4,35, 0.01)
mat.y <- m(y,theta.y.1,theta.y.2)
contour(theta.y.1, theta.y.2,mat.y, bty = "l", las =1, ann=F)
mtext(side = 1, text = expression(theta[y1]), line =2)
mtext(side = 2, text = expression(theta[y2]), line =2)
#side on which side of the plot (1=bottom, 2=left, 3=top, 4=right)
theta.z.1 <- seq (1,6,0.005); theta.z.2 <- seq(5,50, 0.01)
mat.z <- m(z,theta.z.1,theta.z.2)
contour(theta.z.1, theta.z.2,mat.z)


theta.y.1.mar <-rowSums(mat.y)/sum(mat.y)
theta.y.1.draws <- sample(theta.y.1, 1000, prob=theta.y.1.mar)
theta.y.2.draws <- numeric(1000)
for (i in 1:1000){
  theta.y.2.draws[i]<- sample(theta.y.2, 1, prob = mat.y[which(theta.y.1 == theta.y.1.draws[i]),])
}
theta.y.1.draws <- jitter(theta.y.1.draws, amount = 0.005/2)
theta.y.2.draws <- jitter(theta.y.2.draws, amount = 0.01/2)

theta.z.1.mar <-rowSums(mat.z)/sum(mat.z)
theta.z.1.draws <- sample(theta.z.1, 1000, prob=theta.z.1.mar)
theta.z.2.draws <- numeric(1000)
for (i in 1:1000){
  theta.z.2.draws[i]<- sample(theta.z.2, 1, prob = mat.z[which(theta.z.1 == theta.z.1.draws[i]),])
}
theta.z.1.draws <- jitter(theta.z.1.draws, amount = 0.005/2)
theta.z.2.draws <- jitter(theta.z.2.draws, amount = 0.01/2)


## mu.y - mu.z
mu.y <- theta.y.1.draws/(theta.y.2.draws+theta.y.1.draws)
mu.z <- theta.z.1.draws/(theta.z.2.draws+ theta.z.1.draws)
tikz('histdiff.tex', standAlone=TRUE, width=4, height=4)
hist(mu.y- mu.z, breaks =100, prob =T, yaxt = "n", ylab = NULL, main = NULL, ann =F)
#mtext(text =expression(mu[y]-mu[z]), side =1 , line =2)
mtext(side=1, text="$\\mu_y-\\mu_z$", line=2)
dev.off()
tools::texi2dvi('histdiff.tex', pdf=TRUE)

tikz('histdiff`.tex', standAlone=TRUE, width=4, height=4)
hist(rnorm(10000), breaks =100, prob =T, yaxt = "n", 
     ylab = NULL, main = NULL, ann =F)
#mtext(text =expression(mu[y]-mu[z]), side =1 , line =2)
mtext(side=1, text="$\\mu_y-\\mu_z$", line=2)
dev.off()
tools::texi2dvi('histdiff`.tex', pdf=TRUE)

