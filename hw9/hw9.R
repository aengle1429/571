m = 10^5
sigmoid <- function(x)
{
  return(1/(1+exp(-x)))
}
logit <- function(p)
{
  return( log(p/(1-p))
          )
}
vals = matrix(0,nrow=m,ncol=2)
vals2 = vals
x = c()
y = c()
for(i in 1:m)
{
  beta1 = rnorm(1)
  temp = rbinom(1,1,0.5)
  if(temp == 0)
    bi = 0
  else bi = -beta1
  bi2 = rnorm(1)
  x = append(x,beta1)
  y = append(y,bi)
  vals[i,1] = sigmoid(bi)
  vals[i,2] = sigmoid(bi+beta1)
  vals2[i,1] = sigmoid(bi2)
  vals2[i,2] = sigmoid(bi2+beta1)
}

joint <- function(b,beta1)
{
  if(b == 0 || b == -beta1)
  {
    ans = 0.5*dnorm(beta1)
  }
  else ans = 0
  ans
}

Jacobian <- function(mu0, mu1)
{
  return(
    1/(mu0-mu0^2) * 1/(mu1-mu1^2)
  )
}

mesh = 0.001
mu0 <- mu1 <- seq(0+mesh,1-mesh,by=mesh)
z <- matrix(0,nrow=length(mu0),ncol=length(mu0))
w <- z
for(i in 1:length(mu0))
{
  for(j in 1:length(mu0))
  {
    z[i,j] = joint(logit(mu0[i]),logit(mu1[j])-logit(mu0[i]))*Jacobian(mu0[i],mu1[j])
    w[i,j] = dnorm(logit(mu0[i]))*dnorm(logit(mu1[j])-logit(mu0[i]))*Jacobian(mu0[i],mu1[j])
  }
}


contour(mu0,mu1,z,main="Contours from Joint Denisity", xlab="mui0", ylab="mui1")
plot(vals[,1],vals[,2], main="Realizations of mui0, mui1 from 10^5 data points")
contour(mu0,mu1,w,main="Contours from Joint Denisity", xlab="mui0", ylab="mui1")
qplot(vals2[,1],vals2[,2], main="Realizations of mui0, mui1 from 10^5 data points",alpha=I(1/10),xlab="mui0", ylab="mui1")

setwd("~/Dropbox/UW2015-2016/Win2016/571/hw9")
library(lme4)
dat1 <- read.table("pairsdata1.txt", header=T)
dat2 <- read.table("pairsdata2.txt", header=T)
m1 <- glmer(y~x+(1|id), data=dat1, family=binomial)
m2 <- glmer(y~x+(1|id), data=dat2, family=binomial)

pair1 <- dat1[order(dat1$id,dat1$x),]
pair2 <- dat2[order(dat2$id,dat2$x),]

n01 <- sum(pair1$y[seq(1, nrow(pair1), 2)] == 0 & pair1$y[seq(2, nrow(pair1), 2)] == 1)
n10 <- sum(pair1$y[seq(1, nrow(pair1), 2)] == 1 & pair1$y[seq(2, nrow(pair1), 2)] == 0)
n00 <- sum(pair1$y[seq(1, nrow(pair1), 2)] == 0 & pair1$y[seq(2, nrow(pair1), 2)] == 0)
n11 <- sum(pair1$y[seq(1, nrow(pair1), 2)] == 1 & pair1$y[seq(2, nrow(pair1), 2)] == 1)

#first
cbeta1 <- log(n01/n10) #1.306252

n01 <- sum(pair2$y[seq(1, nrow(pair2), 2)] == 0 & pair2$y[seq(2, nrow(pair2), 2)] == 1)
n10 <- sum(pair2$y[seq(1, nrow(pair2), 2)] == 1 & pair2$y[seq(2, nrow(pair2), 2)] == 0)
n00 <- sum(pair2$y[seq(1, nrow(pair2), 2)] == 0 & pair2$y[seq(2, nrow(pair2), 2)] == 0)
n11 <- sum(pair2$y[seq(1, nrow(pair2), 2)] == 1 & pair2$y[seq(2, nrow(pair2), 2)] == 1)

cbeta1 <- log(n01/n10) #1.306252 both times
