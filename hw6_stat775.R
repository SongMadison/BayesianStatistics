## after building a bayesian model, it is very important to test the sensitivity of the model
## and figure the potential limitation of given model in interpreting the data.

# this scirpts means to calculate

#there are total four models:
## check : (1) independent poisson distribution, (2) there is no time trend

set.seed(100)
t <- 1976:1985
acc <- c(24,25,31,31,22,21,26,20,16,22)
y.cor <- cor(acc[-1],acc[-10])
y.p <- summary(lm(acc~t))$coef[8]

theta <- rgamma(1000,sum(acc)+1, rate =10)
cor.rep <- numeric(1000)
p.rep <- numeric(1000)
for( i in 1:1000){
    y <- rpois(10, theta[i])
    cor.rep[i]<- cor(y[-1],y[-10])
    p.rep[i] <- summary(lm(y~t))$coef[8]
}
mean(y.cor<= cor.rep)
hist(cor.rep, breaks =50, prob = T, xlim= c(-1,1), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.cor)
mtext(side =1, text = "Replicated Correlations", line =2)

mean(y.p >= p.rep)
hist(p.rep, breaks =50, prob = T, xlim= c(0,1), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.p)
mtext(side =1, text = "Replicated $p$-values", line =2)


## model 2
death <- c(734,516,754,877,814,362,764,809,223,1066)
death.rate <- c(.19,.12,.15,.16,.14,.06,.13,.13,.03,.15)
x <- death/death.rate*10^8


y.cor <- cor(acc[-1],acc[-10])
y.p <- summary(lm(acc~t))$coef[8]
theta <- rgamma(1000, 238, rate =sum(x))
for ( i in 1:1000){
    y <- rpois(10, theta[i]*x)
    cor.rep[i]<- cor(y[-1],y[-10])
    p.rep[i] <- summary(lm(y~t))$coef[8]
}

mean(y.cor<= cor.rep)
hist(cor.rep, breaks =50, prob = T, xlim= c(-1,1), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.cor)
mtext(side =1, text = "Replicated Correlations", line =2)

mean(y.p >= p.rep)
hist(p.rep, breaks =50, prob = T, xlim= c(0,1), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.p)
mtext(side =1, text = "Replicated $p$-values", line =2)




### consider death as response
#model 3

y.cor <- cor(death[-1],death[-10])
y.p <- summary(lm(death~t))$coef[8]

theta <- rgamma(1000,sum(death), rate =10)
cor.rep <- numeric(1000)
p.rep <- numeric(1000)
for( i in 1:1000){
    y <- rpois(10, theta[i])
    cor.rep[i]<- cor(y[-1],y[-10])
    p.rep[i] <- summary(lm(y~t))$coef[8]
}
mean(y.cor<= cor.rep)
hist(cor.rep, breaks =50, prob = T, xlim= c(-1,1), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.cor)
mtext(side =1, text = "Replicated Correlations", line =2)

mean(y.p >= p.rep)
hist(p.rep, breaks =50, prob = T, xlim= c(0,1), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.p)
mtext(side =1, text = "Replicated $p$-values", line =2)



## model 4
theta <- rgamma(1000,sum(death), rate =sum(x))
cor.rep <- numeric(1000)
p.rep <- numeric(1000)
for( i in 1:1000){
    y <- rpois(10, theta[i]*x)
    cor.rep[i]<- cor(y[-1],y[-10])
    p.rep[i] <- summary(lm(y~t))$coef[8]
}
mean(y.cor<= cor.rep)
hist(cor.rep, breaks =50, prob = T, xlim= c(-1,1), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.cor)
mtext(side =1, text = "Replicated Correlations", line =2)

mean(y.p >= p.rep)
hist(p.rep, breaks =50, prob = T, xlim= c(0,0.003), yaxt ="n", ylab = NULL, main = NULL, ann =F)
abline(v=y.p)
mtext(side =1, text = "Replicated $p$-values", line =2)
