#ex2 redo the experiment about model selection using AIC, BIC
      # more understanding on AIC and BIC
data(mtcars)
str(mtcars)

model1 <- lm (mpg ~ cyl+hp+gear, data =mtcars)
AIC(model1)
sum(model1$residuals^2)




#ex 13 chapter 5
set.seed(10)
mar.post.unnorm <- function(u,v,y,n)
    prod(beta(exp(u+v)/(1+exp(u))+y, exp(v)/(1+exp(u))+n-y) / beta(exp(u+v)/(1+exp(u)), exp(v)/(1+exp(u)))) *exp(u-v/2)/(1+exp(u))^2

m <- function(a,b){
    mat <- matrix(nrow =length(a),ncol=length(b))
    for (i in 1:length(a)){
        for (j in 1:length(b)){
            mat[i,j]<- mar.post.unnorm(a[i],b[j],y,n)
        }
    }
    return(mat)
}

y <- c(16,9,10,13,19,20,18,17,35,55)
n <- c(74,99,58,70,122,77,104,129,308,119)

u <- seq (-1.8, -0.8, 0.001); v <- seq (1.5,4,0.001)
mat <- m(u,v)
mat <- mat/(sum(mat)*0.001*0.001)   # normlizing
contour(u,v,mat, bty = "l",las =1, ann =F)
mtext(side =1, text= expression(log(alpha/beta)),line =2)
mtext(side =2, text= expression(log(alpha+beta)),line =2)


#latex graphics for R
if (!require("tikzDevice")){   
    install.packages("tikzDevice") 
    require("tikzDevice")
}

tikz('histdiff.tex', standAlone=TRUE, width=4, height=4)
contour(u,v,mat, bty = "l",las =1, ann =F)
mtext(side =1, text= "$\\log(\\alpha/\\beta)$",line =2)
dev.off()
tools::texi2dvi('histdiff.tex', pdf=TRUE)


u.mar <- rowSums(mat)/sum(mat)
u.draws <- sample(u,1000,prob = u.mar)
v.draws <- numeric(1000)
for ( i in 1:1000){
    v.draws[i] <- sample(v,1,prob = mat[which(u==u.draws[i]),])
}
u.draws <- jitter(u.draws, amount =0.001/2)
v.draws <- jitter(v.draws, amount =0.001/2)
plot(u.draws,v.draws,pch =20, las =1,bty = "l",ann=F)



#draw smples from joint posterior distribution of parameters and hyperparamenters
alpha.draws <- numeric(1000)
beta.draws <- numeric(1000)
theta.draws <- matrix(data= NA, nrow =1000, ncol =10)


for ( i in 1:1000){
    alpha.draws[i]<- exp(u.draws[i]+v.draws[i])/(1+exp(u.draws[i]))
    beta.draws[i] <- exp(v.draws[i])/(1+exp(u.draws[i]))
    for (j in 1:10){
        theta.draws[i,j]<- rbeta(1,alpha.draws[i]+y[j], beta.draws[i]+n[j]-y[j])  
    }
}

theta.med <- apply(X = theta.draws, MARGIN =2, median)
theta.limits <- apply(X= theta.draws,MARGIN =2, FUN = function (x) quantile(x,probs=c(0.025,0.975)))

plot(y/n, theta.med, pch =20, las = 1, bty="l", ann =F, xlim =c(0,0.5), ylim= c(0,0.5))
abline(0,1)
for( i in 1:10){
    lines(rep(y[i]/n[i], 2),theta.limits[,i],lty=3)
}
legend("topleft", bty="n", pch =20, legend="Median")

## 95% posterior confidence interval for the average underlying proportion of trafic that is bicycle.
theta.means <- apply(X=theta.draws, MARGIN =1, mean)
( theta.conf <- quantile(theta.means,probs = c(0.025,0.975)) )

theta.pop.means <- alpha.draws/(alpha.draws+beta.draws)
quantile(theta.pop.means,probs = c(0.025,0.975))

##e 
y.post <- numeric(1000)
for( i in 1:1000){
    y.post[i] <- rbinom(n = 1,size = 100,prob = theta.means[i])
}
(quantile(y.post, probs =c(0.025, 0.975) ))

y.post <- numeric(1000)
for( i in 1:1000){
    y.post[i] <- rbinom(n = 1,size = 100,prob = theta.pop.means[i])
}
quantile(y.post, probs =c(0.025, 0.975) )

# correct way.
y.post <- numeric(1000)
for( i in 1:1000){
    y.post[i] <- rbinom(n = 1,size = 100,prob = rbeta(1, alpha.draws[i],beta.draws[i]))
}
quantile(y.post, probs =c(0.025, 0.975) )

##

