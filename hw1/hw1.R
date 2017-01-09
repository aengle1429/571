library(sandwich)
setwd("~/Dropbox/UW2015-2016/Win2016/571/hw1")
##
##problem 1
##
dat1 <- read.table("rats.csv", header=TRUE, sep=',')
attach(dat1)
m1 <- glm(y~x, family=Gamma(link="log")) #we have a log link log(1/lambda_i) = beta_0 + beta_1x_i
sand <- vcovHC(m1,type="HC0")
beta <- m1$coefficients
n = length(x)
X = cbind(rep(1,n),x)
temp = exp(X%*%(-beta))
A = (-1/n)*matrix(c(y%*%temp, (x*y)%*%temp, (x*y)%*%temp, sum(x^2*y*temp)
  ), nrow=2, ncol=2
)

fun <- function(x,y)
{
  return(
    matrix(c(
      -1+y*exp(-beta[1]-x*beta[2]), -x+x*y*exp(-beta[1]-x*beta[2])
    ))
  )
}

B = matrix(0,nrow=2,ncol=2)
for(i in 1:n)
{
  B = B + fun(x[i],y[i])%*%t(fun(x[i],y[i]))
}
B = B/n

yum = solve(A)%*%B%*%solve(A)/n
yum2 = solve(Ah)%*%B%*%solve(Ah)/n
detach(dat1)

##
##problem 3
##
library(MASS)
n = dim(cars)[1]
R <- function(rho)
{
  temp = matrix(0, nrow=n, ncol=n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      temp[i,j] = rho^(abs(i-j))
    }
  }
  return(temp)
}

#|rho| < 1 for full column rank
X = cbind(rep(1,n), cars$dist)
probs <- function(sigmasq,rho){
  countint = 0
  countslope = 0
  for(i in 1:10000)
  {
    y = X%*%beta + mvrnorm(1,rep(0,n),sigmasq*R(rho))
    m1 <- lm(y~cars$dist)
    temp = confint.default(m1)
    if(
      findInterval(beta[1], temp[1,])==1 #is beta1 in the 95% interval?
      )
    {
      countint = countint + 1
    }
    if(
      findInterval(beta[2], temp[2,])==1 #is beta2 in the 95% interval?
    )
    {
      countslope = countslope + 1
    }
  }
  return(c(countint/10000
  ,countslope/10000))
}

set.seed(1)
sigmas <- c(0.1,1,10,100)
rhos <- c(-0.9,-0.5,-0.25,0,0.25,0.5,0.9)
beta = c(1,-1)
for(sig in sigmas){
  for(rho in rhos)
  {
    cat("coverage for ", sig, rho, probs(sig,rho),"\n")
  }
}  

beta = c(5,30)
for(sig in sigmas){
  for(rho in rhos)
  {
    cat("coverage for ", sig, rho, probs(sig,rho),"\n")
  }
}  


probs(1,0.5)


H = X%*%solve(t(X)%*%X)%*%t(X)
sum(diag((diag(n)-H)%*%R(0.9))) #

##
##p4
##
p4 <- read.table("evggfpd30.csv", header=TRUE, sep=",")
library(ggplot2)
m <- dim(p4)[1] #this many people
n <- dim(p4)[2]

prop = rep(0,n-1)
for(i in 1:31)
{
  prop[i] = m - sum(p4[,i]==p4[,i+1])
}

prop = prop / m
dftemp = data.frame(cbind(1:31,prop))
plotp1 <- ggplot(dftemp, aes(x=V1, y=prop))
plotp1 + geom_point() + geom_line() + ylim(c(0,1)) + xlab("Time") + ylab("Proportion") + ggtitle("Proportion of People Moving Between Pairs of Health States")

##p4p2

countsPerYear <- apply(p4,2,table) #compute counts of levels per year

zeros = as.matrix(p4) == matrix(0,nrow=m,ncol=n)
ones = as.matrix(p4) == matrix(1,nrow=m,ncol=n)
twos = as.matrix(p4) == matrix(2,nrow=m,ncol=n)
threes = as.matrix(p4) == matrix(3,nrow=m,ncol=n)
fours = as.matrix(p4) == matrix(4,nrow=m,ncol=n)
fives = as.matrix(p4) == matrix(5,nrow=m,ncol=n)
dead = cbind(unname(apply(zeros,2,sum)),rep(0,32))
poor = cbind(unname(apply(ones,2,sum)),rep(1,32))
fair = cbind(unname(apply(twos,2,sum)),rep(2,32))
good = cbind(unname(apply(threes,2,sum)),rep(3,32))
vgood = cbind(unname(apply(fours,2,sum)),rep(4,32))
eggs = cbind(unname(apply(fives,2,sum)),rep(5,32))


p4p2 <- data.frame(cbind(rep(1:32,6),rbind(dead, poor, fair, good, vgood, eggs)), row.names=NULL)
plot <- ggplot(p4p2, aes(X1, X2, color=factor(X3))) + geom_point() + geom_line() + xlab("Time (years)") + ylab("Absolute Count")
plot <- plot + ggtitle("Absolute Number of People in Different Health States Over Time")


##p4p3
set.seed(1)
people <- sample(seq(1,m), 100, replace=TRUE)
temp = as.matrix(p4[people,])
vec = c()
vec2 = c()
for(i in 1:100)
{
  vec = append(vec,temp[i,])
  vec2 = append(vec2, rep(i,32))
}
vec0 = rep(1:32,100)
df <- data.frame(cbind(vec0,vec,vec2), row.names=NULL)


ggplot(df, aes(vec0,vec)) + geom_point() + geom_line() + facet_wrap(~vec2)