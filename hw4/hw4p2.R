library(geeM)
B = 10^4
n = 50
sigma = 1
tau = 1
beta0 = 0
beta1 = 0.8
set.seed(100)
count = 0
ests <- matrix(0,nrow=B,ncol=2)
for(v in 1:B)
{
  nis <- rpois(n,3) + 1
  Z <- matrix(0,nrow=n, ncol=max(nis))
  X <- matrix(0,nrow=n, ncol=max(nis))
  for(i in 1:n)
  {
    for(j in 1:nis[i])
    {
      Z[i,j] = rgamma(1,shape=3,rate=2)
    }
  }
  Z[Z == 0] <- NA
  for(i in 1:n)
  {
    X[i,1:nis[i]] = cumsum(Z[i,1:nis[i]])
  }
  X[X == 0] <- NA
  
  b = rnorm(n,mean=0,sd=tau)
  mu = b + beta0 + beta1*X
  Y <- matrix(0,nrow=n, ncol=max(nis))
  for(i in 1:n)
  {
    for(j in 1:nis[i])
    {
      Y[i,j] = rnorm(1,mean=mu[i,j],sd=sigma)
    }
  }
  Y[Y==0] <- NA
  
  ys = as.vector(t(Y))
  xs = as.vector(t(X))
  simulation = data.frame(ys[!is.na(ys)], xs[!is.na(xs)], rep(1:n,times=nis))
  names(simulation) <- c("y", "x", "id")
  m1 <- geem(y~x, id=id, data=simulation, corstr="independence")
  ci <- 1.96*summary(m1)$se.robust[2]*c(-1,1) + m1$beta[2]
  ests[v,1] = m1$beta[2]
  ests[v,2]=summary(m1)$se.robust[2]
  if(findInterval(beta1,ci) == 1)
  {
    count = count + 1
  }
}
hist(ests[,1],main="Histogram of Slope Estimates, n=50, B=10^4", xlab=expression(hat(beta)[50]))
abline(v=beta1, col="red")
legend("topright",legend=expression(paste("True ", beta)),lty=1,lwd=1,col="red")

estStdErr <- mean(ests[,2])


hist(ests[,2],xlab="Estimated Std Err", main="Histogram of Estimated Standard Errors, n=50, B=10^4")
abline(v=estStdErr, col="red")
legend("topright",legend=expression(paste("Est. of True Std Err")),lty=1,lwd=1,col="red")

plot(qnorm(ppoints(B)),quantile(ests[,1],probs=ppoints(B)),
     main=expression(paste("QQ plot of B=10^4 Replicates of ", hat(beta)[50])),
     xlab="Theoretical Quantiles", ylab="Sample Quantiles")
abline(a=beta1, b=estStdErr,col="red",lwd=2)
legend("bottomright", "Intercept = True Parameter\n Slope=Est of True Std Err",lty=1,lwd=2,col="red")

#plot for joint density
bivariate = data.frame(c(ests[,1],beta1), c(ests[,2],estStdErr),as.factor(rep(c(1,2),times=c(B,1))))
names(bivariate) <- c("beta", "stderr", "id")
ggplot(bivariate, aes(beta, stderr, color=id)) + geom_point(alpha=1/2) + ggtitle("Estimates from B=10^4 Replications")

plot(ests[,1],ests[,2],col="blue")
points(c(beta1,estStdErr),col="red")

# > count
# [1] 9129