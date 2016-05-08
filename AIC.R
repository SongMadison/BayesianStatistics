data(mtcars)
str(mtcars)

model1 <- lm (mpg ~ cyl+hp+gear, data =mtcars)
AIC(model1)


n = nrow(mtcars)
RSS = sum(model1$residuals^2)
n * log(RSS/n) + 2*4

n/2 * log(RSS/n) + 2*4 -n/2*log(2*pi) - n/2*log ()