library(ggplot2)
library(gee)
library(geeM)
n = 10000
beta1 = 5
beta2 = seq(-5,5,0.5)
#something in [-5,5]
beta1s <- matrix(0,nrow=1,ncol=length(beta2))
beta1sAR1 <- matrix(0,nrow=1,ncol=length(beta2))

set.seed(1000)
beta2=5
for(k in 1:length(beta2))
{
  X = matrix(rep(1:4,n),nrow=n, ncol=4, byrow=TRUE)
  mu = beta1 * X + beta2[k]*(X==2)-beta2[k]*(X==3)
  Y = matrix(0,nrow=n, ncol=4)
  for(i in 1:n)
  {
    for(j in 1:4)
    {
      Y[i,j] = rnorm(1,mean=mu[i,j], sd=1)
    }
  }
  ys = as.vector(t(Y))
  xs = as.vector(t(X))
  simulations = data.frame(ys,xs,rep(1:n,each=4))
  names(simulations) <- c("y", "x", "id")
  m1 <- geem(y~x, id=id, data=simulations, corstr = "independence")
  beta1s[1,k] = m1$beta[2]
}

#now do it for working cor = AR-1
beta2=5
for(k in 1:length(beta2))
{
  X = matrix(rep(1:4,n),nrow=n, ncol=4, byrow=TRUE)
  mu = beta1 * X + beta2[k]*(X==2)-beta2[k]*(X==3)
  Y = matrix(0,nrow=n, ncol=4)
  for(i in 1:n)
  {
    for(j in 1:4)
    {
      Y[i,j] = rnorm(1,mean=mu[i,j], sd=1)
    }
  }
  ys = as.vector(t(Y))
  xs = as.vector(t(X))
  simulations = data.frame(ys,xs,rep(1:n,each=4))
  names(simulations) <- c("y", "x", "id")
  m1 <- geem(y~x, id=id, data=simulations, corstr = "ar1")
  beta1sAR1[1,k] = m1$beta[2]
}
betaBox <- data.frame(c(beta1s,beta1sAR1),rep(seq(-5,5,0.5),2),as.factor(rep(c("Indp","AR1"),each=length(beta2))))
names(betaBox) <- c("slopes", "true", "id")
ggplot(betaBox, aes(x=true, y=slopes,col=id))  + geom_line() + ggtitle(expression(paste("Slopes against ", beta[2], " Values")))+
  xlab("Choice of Beta 2") + ylab("Limiting Value of Slope")
