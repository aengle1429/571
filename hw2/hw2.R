library(rootSolve)
library(ggplot2)
setwd("~/Dropbox/UW2015-2016/Win2016/571/hw2")
setwd("C:/Users/aengl_000/Dropbox/UW2015-2016/Win2016/571/hw2")
lep <- read.table("leprosy.txt", header=TRUE, sep=" ")
#treatment 1 - antibiotic A
#treatment 2 - antibiotic B
#treatment 3 - placebo
#count1: Pre-treatment count of bacilli, at six sites of the body where the bacilli tend to congregate
#count2: Post-treatment count of bacilli (lower is better)
#severe: Indicator of severe disease, prior to the trial (0=No, 1=Yes)
attach(lep)
n = dim(lep)[1]
ni = 2 #same for all i
p = 4 #four parameters

X = matrix(0,nrow=60,ncol=4) #sum(ni) x p is 60 x 4 stacked matrix
X[seq(1,59,by=2),] <- rep(c(1,0,0,0),each=30)
X[seq(2,60,by=2),] <- rep(c(0,1,0,0),each=30)
for(i in 1:30)
{
  if(trt[i] == 1)
  {
    X[2*i,3] = 1
  }
  if(trt[i] == 2)
  {
    X[2*i,4] = 1
  }
}
#Y = rbind(count1,count2)
Y = rbind(log(count1)),log(count2))
est <- function(beta){
  temp = 0
  count = 1
  for(i in seq(1,59,by=2))
  {
    temp = temp + t(X[c(i,i+1),])%*%(t(t(Y[,count]))-X[c(i,i+1),]%*%beta)
    count = count +1
  }
  return(
    (1/n)*temp
  )
}
beta = multiroot(est,start=c(10,0,0,0))$root
theta = c(beta[1],beta[2]/beta[1] ,(beta[2]+beta[3])/beta[2] ,(beta[2]+beta[4])/beta[2])

set.seed(1)
boot <- function(B)
{
  betas <- matrix(0, nrow=B, ncol=4)
  for(b in 1:B)
  {
    ind <- sample(1:30, 30, replace=TRUE)
    tempdat <- lep[ind,]
    X = matrix(0,nrow=60,ncol=4) #sum(ni) x p is 60 x 4 stacked matrix
    X[seq(1,59,by=2),] <- rep(c(1,0,0,0),each=30)
    X[seq(2,60,by=2),] <- rep(c(0,1,0,0),each=30)
    for(i in 1:30)
    {
      if(tempdat$trt[i] == 1)
      {
        X[2*i,3] = 1
      }
      if(tempdat$trt[i] == 2)
      {
        X[2*i,4] = 1
      }
    }
    Y = rbind(tempdat$count1,tempdat$count2)
    est <- function(beta){
      temp = 0
      count = 1
      for(i in seq(1,59,by=2))
      {
        temp = temp + t(X[c(i,i+1),])%*%(t(t(Y[,count]))-X[c(i,i+1),]%*%beta)
        count = count +1
      }
      return(
        (1/n)*temp
      )
    }
    beta = multiroot(est,start=c(0,0,0,0))$root
    betas[b,] <- beta
    #thetas[b,] <- c(beta[1],beta[2]/beta[1] ,(beta[2]+beta[3])/beta[2] ,(beta[2]+beta[4])/beta[2])
  }
  return(betas)
}
betas <- boot(10000)
dat = data.frame(bs=c(betas[,1], betas[,2], betas[,3], betas[,4]), lab=as.factor(rep(1:4, each=10^4)))

ggplot(dat, aes(x=bs, fill=lab)) +
  geom_histogram(alpha=0.2, position="identity", binwidth=0.5) +
  xlab("Sample Values")+
  ggtitle("Marginal Histograms of Betas")

theta1 <- betas[,1]
theta2 <- betas[,2]/betas[,1]
theta3 <- (betas[,2] + betas[,3])/betas[,2]
theta4 <- (betas[,2] + betas[,4])/betas[,2]

dat2 = data.frame(bs=c(theta2,theta3,theta4), theta=as.factor(rep(2:4, each=10^4)))

ggplot(dat2, aes(x=bs, fill=theta)) +
  geom_histogram(binwidth=0.1) + facet_grid(theta~.)+
  xlab("Sample Values")+
  ggtitle("Marginal Histograms of Theta 1,2,3")


dat3 = data.frame(bs=c(log(theta3),log(theta4)), theta=as.factor(rep(3:4, each=10^4)))
ggplot(dat3, aes(x=bs, fill=theta)) +
  geom_histogram(binwidth=0.1) + facet_grid(theta~.)+
  xlab("Sample Values")+
  ggtitle("Marginal Histograms of Log-transformed Theta 2,3")

qqnorm(log(theta3))
qqline(log(theta3))

qqnorm(log(theta4))
qqline(log(theta4))


dat4 = data.frame(bs=c(sqrt(theta3),sqrt(theta4)), theta=as.factor(rep(3:4, each=10^4)))
ggplot(dat3, aes(x=bs, fill=theta)) +
  geom_histogram(binwidth=0.1) + facet_grid(theta~.)+
  xlab("Sample Values")+
  ggtitle("Marginal Histograms of SquareRoot-transformed Theta 2,3")
qqnorm(sqrt(theta3))
qqline(sqrt(theta3))

qqnorm(sqrt(theta4))
qqline(sqrt(theta4))

CI1 <- quantile(theta1,probs=c(0.025,0.5,0.975))
CI2 <- quantile(theta2,probs=c(0.025,0.5, 0.975))
CI3 <- quantile(sqrt(theta3),probs=c(0.025,0.5,0.975))
CI4 <- quantile(sqrt(theta4),probs=c(0.025,0.5,0.975))
#p2p4
J3 <- function(x)
{
  p = length(x)
  return(
    matrix(c(1,0,0,0,
             -x[2]/x[1]^2,1/x[1],0,0,
             0,-x[3]/(2*x[2]^2*sqrt((x[2]+x[3])/x[2])),1/(2*x[2]*sqrt((x[2]+x[3])/x[2])),0,
             0,-x[4]/(2*x[2]^2*sqrt((x[2]+x[4])/x[2])),0,1/(2*x[2]*sqrt((x[2]+x[4])/x[2]))), nrow=p, ncol=p, byrow=TRUE)
  )
}
Bhat = 0
count = 1
for(i in seq(1,59,by=2))
{
  Bhat = Bhat + t(X[c(i,i+1),])%*%(t(t(Y[,count]))-X[c(i,i+1),]%*%beta)%*%t(t(X[c(i,i+1),])%*%(t(t(Y[,count]))-X[c(i,i+1),]%*%beta))
  count = count +1
}
Bhat = Bhat/n

Ahat = 0
for(i in seq(1,59,by=2))
{
  Ahat = Ahat - t(X[c(i,i+1),])%*%X[c(i,i+1),]
}
Ahat = Ahat/n

thetanew2 = c(beta[1],beta[2]/beta[1] ,sqrt((beta[2]+beta[3])/beta[2]) ,sqrt((beta[2]+beta[4])/beta[2]))
sandNEW2 = J3(beta)%*%solve(Ahat)%*%Bhat%*%solve(t(Ahat))%*%t(J3(beta))/n
CIsand1NEW2 <- thetanew2[1] + 1.96*sqrt(sandNEW2[1,1])*c(-1,1)
CIsand2NEW2 <- thetanew2[2] + 1.96*sqrt(sandNEW2[2,2])*c(-1,1)
CIsand3NEW2 <- thetanew2[3] + 1.96*sqrt(sandNEW2[3,3])*c(-1,1)
CIsand4NEW2 <- thetanew2[4] + 1.96*sqrt(sandNEW2[4,4])*c(-1,1)


#not very good...
theta <- c(beta[1],beta[2]/beta[1] ,(beta[2]+beta[3])/beta[2] ,(beta[2]+beta[4])/beta[2])
sand = J(beta)%*%solve(Ahat)%*%Bhat%*%solve(t(Ahat))%*%t(J(beta))/n
CIsand1 <- theta[1] + 1.96*sqrt(sand[1,1])*c(-1,1)
CIsand2 <- theta[2] + 1.96*sqrt(sand[2,2])*c(-1,1)
CIsand3 <- theta[3] + 1.96*sqrt(sand[3,3])*c(-1,1)
CIsand4 <- theta[4] + 1.96*sqrt(sand[4,4])*c(-1,1)

sandt <- solve(Ahat)%*%Bhat%*%solve(t(Ahat))/n
sqrt(diag(sandt))

#these below are from untransformed, very good
CIsandt1 <- beta[1] + 1.96*sqrt(sandt[1,1])*c(-1,1)
CIsandt2 <- beta[2] + 1.96*sqrt(sandt[2,2])*c(-1,1)
CIsandt3 <- beta[3] + 1.96*sqrt(sandt[3,3])*c(-1,1)
CIsandt4 <- beta[4] + 1.96*sqrt(sandt[4,4])*c(-1,1)

CIt1 <- quantile(betas[,1],probs=c(0.025,0.975))
CIt2 <- quantile(betas[,2],probs=c(0.025,0.975))
CIt3 <- quantile(betas[,3],probs=c(0.025,0.975))
CIt4 <- quantile(betas[,4],probs=c(0.025,0.975))

# ##new transform of parameters
# J2 <- function(x)
# {
#   p = length(x)
#   return(
#     matrix(c(1,0,0,0,
#              -x[2]/x[1]^2,1/x[1],0,0,
#              0,-x[3]/(x[2]^2+x[2]*x[3]),1/(x[2]+x[3]),0,
#              0,-x[4]/(x[2]^2+x[2]*x[4]),0,1/(x[2]+x[4])), nrow=p, ncol=p, byrow=TRUE)
#   )
# }
# thetanew = c(beta[1],beta[2]/beta[1] ,log((beta[2]+beta[3])/beta[2]) ,log((beta[2]+beta[4])/beta[2]))
# sandNEW = J2(beta)%*%solve(Ahat)%*%Bhat%*%solve(t(Ahat))%*%t(J2(beta))/n
# CIsand1NEW <- thetanew[1] + 1.96*sqrt(sandNEW[1,1])*c(-1,1)
# CIsand2NEW <- thetanew[2] + 1.96*sqrt(sandNEW[2,2])*c(-1,1)
# CIsand3NEW <- thetanew[3] + 1.96*sqrt(sandNEW[3,3])*c(-1,1)
# CIsand4NEW <- thetanew[4] + 1.96*sqrt(sandNEW[4,4])*c(-1,1)
# 
# theta1NEW <- betas[,1]
# theta2NEW <- betas[,2]/betas[,1]
# theta3NEW <- log((betas[,2] + betas[,3])/betas[,2])
# theta4NEW <- log((betas[,2] + betas[,4])/betas[,2])
# CI1new <- quantile(theta1NEW,probs=c(0.025,0.975))
# CI2new <- quantile(theta2NEW,probs=c(0.025,0.975))
# CI3new <- quantile(theta3NEW,probs=c(0.025,0.975))
# CI4new <- quantile(theta4NEW,probs=c(0.025,0.975))