

library(reshape2)
library(ggplot2)
source("R/utilities_wv.R")
Rcpp::sourceCpp('src/utilities.cpp')
library(MASS)
set.seed(20190531)

X = scale(rnorm(100))*sqrt(100)/sqrt(99); q1 = rep(0, 1000)
for (i in 1:1000){
  Y = mvrnorm(100, rep(0,2), matrix(c(1,0.5, 0.5, 1), nrow=2))
  q1[i] = get_q(Y[,1], Y[,2], X, coef = 0, correction = FALSE)
}
qqplot(q1, rchisq(1000, 1)); abline(0,1,col = 'red')




for (i in 1:1000){
  Y = mvrnorm(100, rep(0,2), matrix(c(2,0.5, 0.5, 1), nrow=2))
  q1[i] = get_q(Y[,1], Y[,2], X, coef = 0, correction = FALSE)
}
qqplot(q1, rchisq(1000, 1)); abline(0,1,col = 'red')



Y = mvrnorm(100, rep(0,2), matrix(c(2,0.5, 0.5, 1), nrow=2))
Y2 = cbind(Y[,1]/sqrt(mean(Y[,1]^2)), Y[,2]/sqrt(mean(Y[,2]^2)))
Y3 = scale(Y) * sqrt(100)/sqrt(99)
plot(Y2[,1] ~ Y3[,1]); abline(0,1,col = 'red')
plot(Y2[,2] ~ Y3[,2]); abline(0,1,col = 'red')
