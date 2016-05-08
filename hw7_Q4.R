


est1 <- numeric(100)
for ( i in 1:100){
    x <- rnorm(100)
    est1[i] <- mean(x>2)
}
hist(est1)
mean(est1)
var(est1)


est2<- numeric(100)
for ( i in 1:100){
    x <- rcauchy(100)
    w <- dnorm(x)/dcauchy(x)
    est2[i] <- sum((x>2)*w)/sum(w)
}
hist(est2)
mean(est2)
var(est2)
