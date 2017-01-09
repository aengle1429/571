setwd("~/Dropbox/UW2015-2016/Win2016/571/hw3")
library(ggplot2)
library(geeM)
library(reshape2)
library(boot)
##prob1
a = qnorm(ppoints(100))
g <- function(x)
{
  return(
    exp(x)/(1+exp(x))
  )
}
n = 10000
X = cbind(rep(1,2),c(0,1))
beta1 = seq(-2,2,0.5)
b1hat = rep(0,length(beta1))
b0hat = b1hat
trueb1hat = b1hat
trueb0hat = b1hat
medianDisc = b1hat
set.seed(1)
for(k in 1:length(beta1))
{
  simulation = matrix(0,nrow=n,ncol=2)
  discpairs = rep(0,n)
  for(i in 1:n){
    N = 100 #total number of matched pairs
    Y = matrix(0,nrow=100,ncol=2)
    m01 = 0 #counts for discordant pairs
    m10 = 0
    for(j in 1:N)
    {
      beta = c(a[j],beta1[k])
      Y[j,1] = rbinom(1,1,g(X%*%beta)[1])
      Y[j,2] = rbinom(1,1,g(X%*%beta)[2])
      if(Y[j,1] == 0 & Y[j,2] == 1)
      {
        m01 = m01 + 1
      }
      if(Y[j,1] == 1 & Y[j,2] == 0)
      {
        m10 = m10 + 1
      }
    }
  discpairs[i] <- m01/m10
  simulation[i,] <- colMeans(Y)
  }
  
  trueb0hat[k] = logit(mean(g(a)))
  trueb1hat[k] = logit(mean(g(a+beta1[k]))) - logit(mean(g(a)))
  
  cat("Beta1 is ", beta1[k], "\n")
  cat("Median of Discordant Pairs is ", median(log(discpairs)), "\n")
  medianDisc[k] = median(log(discpairs))
  cat("First Estimate is ", logit(mean(simulation[,1])), "\n")
  b0hat[k] = logit(mean(simulation[,1]))
  cat("Second Estimate is ", logit(mean(simulation[,2]))-logit(mean(simulation[,1])), "\n")
  b1hat[k] = logit(mean(simulation[,2]))-logit(mean(simulation[,1]))
}

plot(beta1,trueb0hat,ylim=c(-0.01,0.01),col="blue")
points(beta1,b0hat,col="red")

plot(beta1,abs(trueb1hat-b1hat),col="blue")
points(beta1,b1hat,col="red")


hist(simulation[,1],main=expression(paste("Histogram of Estimates ", hat(beta)[0])))
abline(v=mean(g(a)),col="red", lty=4, lwd=3)

hist(simulation[,2], main=expression(paste("Histogram of Case Exposure Proportions")))
abline(v=mean(g(a+beta1)),col="red", lty=4, lwd=3)

hist(logit((simulation[,2]))-logit((simulation[,1])),main=expression(paste("Histogram of Estimates ", hat(beta)[1])))
abline(v=logit(mean(g(a+beta1)))-logit(mean(g(a))),col="red", lty=4, lwd=3)

##########################
###estimates of b0hat#####
#> b0hat
#[1] -0.001900001  0.000344000 -0.002124001  0.001536000  0.001520000  0.000480000
#[7]  0.000428000 -0.001200000 -0.000772000
##########################
##########################



(1-g(a))*(g(a+beta1[1]))
#all probabilities are <0.07
#beta1 = -2
#this is the probability of (0,1) goes into m01 count


(g(a))*(1-g(a+beta1[length(beta1)]))
#all probabilities are <0.07
#beta1 = 2
#this is the probability of (1,0), goes into m10 count


#########
##prob2##
#########
library(nlme)

data(Orthodont, package="nlme")

d4 <- Orthodont # using shorter names
d4$id <- d4$Subject
d4$male <- d4$Sex=="Male"

m1 <- glm(distance~I(age-8)*male,data=d4)
beta = m1$coefficients
betaold = beta
X = model.matrix(m1)

workcor = "AR1"
count = 0
repeat{
  count = count + 1
  phi = (1/(4*27-4))*sum((d4$distance - X%*%beta)^2)
  pearsonResid = phi^(-1/2)*(d4$distance - X%*%beta)
  
  pearson_resid = matrix(pearsonResid,nrow=27,ncol=4,byrow=TRUE)
  
  start = 4*0:26+1
  n = 27
  ni = 4
  p = 4
  
  if(workcor == "AR1"){
  alphaAR1 = 0
    for(i in 1:n){
      for(j in 1:(ni-1)){
        alphaAR1 <- alphaAR1 + pearson_resid[i, j] * pearson_resid[i, j+1]
      }
    }
    alphaAR1 <- 1 / (n*ni - n - p) * alphaAR1
  }
  
  
  if(workcor == "Exchangeable")
  {
  alphaExch = 0
  for(i in 1:n)
  {
    for(j in 1:(ni-1))
    {
      for(k in (j+1):ni)
      {
        alphaExch <- alphaExch + pearson_resid[i, j] * pearson_resid[i, k]
      }
    }
  }
  alphaExch <- 1 / (n / 2 * ni * (ni - 1) - p) * alphaExch
  }
  ##################################
  ##matrices Ri do not vary over i##
  ##################################
  RiExch = matrix(alphaExch,nrow=4,ncol=4)
  diag(RiExch) <- rep(1,4)
  
  RiAR1Inv = matrix(0,nrow=4,ncol=4)
  diag(RiAR1Inv) = 1 + alphaAR1^2
  RiAR1Inv[1,1] = 1
  RiAR1Inv[4,4] = 1
  RiAR1Inv[row(RiAR1Inv)==(col(RiAR1Inv)+1)] = -alphaAR1
  RiAR1Inv[(row(RiAR1Inv)+1)==col(RiAR1Inv)] = -alphaAR1
  RiAR1Inv = RiAR1Inv/(1-alphaAR1^2)
  
  ViInvAR1 = (1/phi)*RiAR1Inv
  ViInvExch = (1/phi)*solve(RiExch)
  
  HessAR1 = matrix(0,nrow=4,ncol=4)
  gradAR1 = matrix(0,nrow=4,ncol=1)
  HessExch = matrix(0,nrow=4,ncol=4)
  gradExch = matrix(0,nrow=4,ncol=1)
  for(i in 1:27)
  {
    Di = X[start[i]:(start[i]+3),]
    
    HessAR1 = HessAR1 + t(Di)%*%ViInvAR1%*%Di
    gradAR1 = gradAR1 + t(Di)%*%ViInvAR1%*%(d4$distance[start[i]:(start[i]+3)] - Di%*%beta)
  
    HessExch = HessExch + t(Di)%*%ViInvExch%*%Di
    gradExch = gradExch + t(Di)%*%ViInvExch%*%(d4$distance[start[i]:(start[i]+3)] - Di%*%beta)
  }
  if(workcor == "Exchangeable"){
  beta = beta + solve(HessExch)%*%gradExch}
  else{
  beta = beta + solve(HessAR1)%*%gradAR1
  }
  if(sum((beta-betaold)^2) < 10^(-6))
  {
    break
  }
  betaold = beta
}

alphaExch - gee2$alpha

alphaAR1 - gee1$alpha

gee1 <- geem(distance~I(age-8)*male, id=id, data=d4, corstr="ar1")
gee2 <- geem(distance~I(age-8)*male, id=id, data=d4, corstr="exchangeable")




#########
##prob3##
#########
lep <- read.table("leprosy.txt", header=TRUE, sep=" ")
attach(lep)
n = dim(lep)[1]
ni = 2 #same for all i
p = 4 #four parameters

lep['severe']<-NULL
lep['subject'] = 1:30
lep_Long = melt(lep, id.vars = c("subject", "trt"))
lep_Long = lep_Long[order(lep_Long$subject),]
lep_Long['pre'] <- rep(c(1,0),30)
lep_Long['post'] <- rep(c(0,1),30)
lep_Long['A'] = as.numeric(lep_Long$trt==1 & lep_Long$post==1)
lep_Long['B'] = as.numeric(lep_Long$trt==2 & lep_Long$post==1)
gee1 <- gee(value~-1+pre+post+A+B, id=subject, data=lep_Long, corstr="independence")
gee2 <- gee(value~-1+pre+post+A+B, id=subject, data=lep_Long, corstr="unstructured")
summary(gee1)
