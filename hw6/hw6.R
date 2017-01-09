library(geeM)
library(gee)
library(reshape2)
library(ggplot2)
library(mitools)
#p1
set.seed(571)
sb = 1
sy = 0.5
beta0 = 1
beta1 = -1
beta2 = 0.5
n = 10^3
betas = matrix(0,nrow=10^3,ncol=3)
count0 = 0
count1 = 0
count2 = 0
for(k in 1:10^3)
{
  b = rnorm(n,mean=0,sd=sb)
  X = rbinom(n,1,0.5)
  y = c()
  for(i in 1:n)
  {
    mu = beta0 + b[i] + beta1*(1:3) + beta2*X[i] 
    y = append(y,rnorm(3, mean = mu, sd = sy))
  }
  dat = data.frame(Y = y,X = rep(X,each=3),t = rep(1:3,n),id = rep(1:n,each=3))
  missing = dat
  inds = 3*(0:(n-1))+1
  for(i in inds)
  {
    for(j in 0:1)
    {
      if(dat$Y[i+j] < 0)
      {
        missing$Y[((i+j+1):(i+2))] <- NA
        break
      }
    }
  }
  model = geem(Y~t+X,id=id,data=dat,corstr="independence")
  model3 = geem(Y~t+X,id=id,data=missing,corstr="independence")
  betas[k,] = model3$beta
  len = 1.96*summary(model3)$se.robust
  if(
    findInterval(beta0,model3$beta[1] + c(-1,1)*len[1]) == 1
  )
  {
    count0 = count0 + 1
  }
  if(
    findInterval(beta1,model3$beta[2] + c(-1,1)*len[2]) == 1
  )
  {
    count1 = count1 + 1
  }
  if(
    findInterval(beta2,model3$beta[3] + c(-1,1)*len[3]) == 1
  )
  {
    count2 = count2 + 1
  }
}
hist(betas[,1],main=expression(paste("Histogram of ", hat(beta)[0])),xlab=expression(hat(beta)[0]))
hist(betas[,2],main=expression(paste("Histogram of ", hat(beta)[1])),xlab=expression(hat(beta)[1]))
hist(betas[,3],main=expression(paste("Histogram of ", hat(beta)[2])),xlab=expression(hat(beta)[2]))
hists = melt(data.frame(cbind(1:10^3,betas)),'X1')
ggplot(hists, aes(x=value)) + geom_histogram(binwidth=0.1) + facet_wrap(~variable)
ggplot(m,aes(value)) + geom_bar(binwidth = 1) + facet_wrap(~variable)
# > count0
# [1] 0
# > count1
# [1] 0
# > count2
# [1] 485
#1b
censor <- function(x,b)
{
  return(
    pnorm((-b-0.5*x)/0.5) + pnorm((-b-0.5*x+1)/0.5) - pnorm((-b-0.5*x)/0.5)*pnorm((-b-0.5*x+1)/0.5)
  )
}
b = seq(-1,1,0.1)
plot(b,1-censor(0,b),main="Conditional Probability of Observing Cluster given X=0",ylab="probability")
plot(b,1-censor(1,b),main="Conditional Probability of Observing Cluster given X=1", ylab="probability")

hist(subset(missing,t==1)$Y)
hist(subset(missing,t==2)$Y)
hist(subset(missing,t==3)$Y)

##1c
n = 10^2
count0 = 0
count1 = 0
count2 = 0
NEWbetas = matrix(0,nrow=1,ncol=3)
for(m in 1)
{
  b = rnorm(n,mean=0,sd=sb)
  X = rbinom(n,1,0.5)
  y = c()
  for(i in 1:n)
  {
    mu = beta0 + b[i] + beta1*(1:3) + beta2*X[i] 
    y = append(y,rnorm(3, mean = mu, sd = sy))
  }
  dat = data.frame(Y = y,X = rep(X,each=3),t = rep(1:3,n),id = rep(1:n,each=3))
  missing = dat
  inds = 3*(0:(n-1))+1
  temp = c()
  for(i in inds)
  {
    for(j in 0:1)
    {
      if(dat$Y[i+j] < 0)
      {
        missing$Y[((i+j+1):(i+2))] <- NA
        break
      }
    }
  }
  model = gee(Y~t+X,id=id,data=dat,corstr="independence")
  model2 = gee(Y~t+X,id=id,data=missing,corstr="independence")
  
  #focus on missing DF
  time1 = subset(missing,t==1)
  miss2 = is.na(missing$Y) & missing$t==2
  miss3 = is.na(missing$Y) & missing$t==3
  clusters = missing[miss2,4] #CLUSTERS MISSING A TIME 2 INDEX
  firstImpY = subset(missing,t==2)$Y
  a = subset(missing,id%in%clusters & t==1) #used for prediction
  mod2 <- lm(firstImpY ~ Y*X, data=time1)
  secondImpY = subset(missing,t==3)$Y
  first2 = subset(missing, t!=3)
  wide = reshape(first2,timevar="t", idvar=c("X","id"),direction="wide")
  mod3 <- lm(secondImpY~Y.1*Y.2*X, data=wide)
  param2 = rnorm(4,mean=summary(mod2)$coefficients[,1],sd=summary(mod2)$coefficients[,2])
  mod2$coefficients = param2
  
  
  param3 = rnorm(8,mean=summary(mod3)$coefficients[,1],sd=summary(mod3)$coefficients[,2])
  mod3$coefficients = param3
  
  sig2 = sqrt(summary(mod2)$sigma^2*rchisq(1,mod2$df.residual)/mod2$df.residual)
  sig3 = sqrt(summary(mod3)$sigma^2*rchisq(1,mod3$df.residual)/mod3$df.residual)
  
  do.one <- function(n)
  {
    miss.i <- missing
  
    mu2 = predict(mod2,a[,-c(3,4)]);
    miss.i[ miss2, "Y"] <- rnorm( sum( miss2), mu2, summary(mod2)$sigma )
    
    
    clusters3 = miss.i[miss3,4]
    b = subset(miss.i,id%in%clusters3 & t!=3)
    blah = reshape(b,timevar="t", idvar=c("X", "id"),direction="wide")
    blah = blah[-2]
    mu3 = predict(mod3,blah);
    
    miss.i[miss3, "Y"] <- rnorm(sum(miss3), mu3, summary(mod3)$sigma)
    
    gee.i <- gee(Y~t+X, id=id,
                  data=miss.i, corstr="independence")
    gee.i
  }
  
  my.mi <- replicate(20, do.one(), simplify=FALSE)
  betas <- MIextract(my.mi, fun=coef)
  vars <- MIextract(my.mi, fun=function(x){as.matrix(x$robust.variance)})
  imputed = MIcombine(betas,vars)
  low = summary(imputed)$`(lower`
  up = summary(imputed)$`upper)`
  NEWbetas[m,] = imputed$coefficients
  if( findInterval(beta0,rbind(low,up)[,1]) == 1)
  {
    count0 = count0 + 1
  }
  if( findInterval(beta1,rbind(low,up)[,2]) == 1)
  {
    count1 = count1 + 1
  }
  if( findInterval(beta2,rbind(low,up)[,3]) == 1)
  {
    count2 = count2 + 1
  }
}

coplot(Y~t|factor(id),data=missing,show.given=F,
       panel=function(x,y,col,...){points(x,y,col=col)
         lines(x,y,col=col)})

matplot(missing$t,missing$Y,type="l",lty=1,
         xlab="time", ylab="Observation", main="Spaghetti Plot of Observations")

mine = NEWbetas[1:1000,]
hist(mine[,1],main=expression(paste("Histogram of ", hat(beta)[0])),xlab=expression(hat(beta)[0]))
hist(mine[,2],main=expression(paste("Histogram of ", hat(beta)[1])),xlab=expression(hat(beta)[1]))
hist(mine[,3],main=expression(paste("Histogram of ", hat(beta)[2])),xlab=expression(hat(beta)[2]))
