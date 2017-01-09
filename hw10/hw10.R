library(geeM)
library(R2WinBUGS)
library(lme4)
library(nlme)
library(mcmc)
library(rjags)
library("INLA")
library("sandwich")
setwd("~/Dropbox/UW2015-2016/Win2016/571/hw10")
setwd("C:\\Users\\aengl_000\\Dropbox\\UW2015-2016\\Win2016\\571\\hw10")
###########
#problem 1#
###########
seiz.l = read.table("seizl.txt", header=TRUE, sep=" ")
Ti <- aggregate( subset(seiz.l, t>1)$y, by=list(seiz.l$id), sum )[,2]
Yi1 <- subset(seiz.l, t==8)$y
trtbas <- subset(seiz.l, t==8)$trt # treatment at time j=1

glm1 <- glm(cbind(Yi1, Ti-Yi1) ~ trtbas , family=binomial)
round( cbind( est=coef(glm1), rob.se=sqrt(diag(vcovHC(glm1, "HC0"))) ) , 2)

Ylong <- rep( rep(c(1,0),each=59), c(Yi1, Ti-Yi1) )
trtbaslong <- rep( c(trtbas,trtbas), c(Yi1, Ti-Yi1) )
glm2 <- glm(Ylong ~ trtbaslong, family=binomial)
round( cbind( est=coef(glm2), rob.se=sqrt(diag(vcovHC(glm2, "HC0"))) ) , 2)

###########
#problem 2#
###########
prob2DF <- read.table("hw10q2.csv", header=TRUE, sep=",")
prob2DF$index = 1:nrow(prob2DF)
lmer1 <- lmer(y~trt+time+trt*time + (1 |id) ,data=prob2DF)
hyperprior <- list(theta=list(prior="loggamma",param=c(0.1,0.1)))
fixedprior <- list( mean.intercept=0, prec.intercept=0.0001,
                    mean=0, prec=0.0001 )
formula1 <- y ~ trt + time + trt*time +
  f(id,model="iid",hyper=hyperprior) #+ f(index,model="iid", hyper=hyperprior)
set.seed(10)
inla1 <- inla(formula1, family="gaussian", data=na.omit(prob2DF), control.fixed=fixedprior, control.family = list(prior="loggamma",param=c(0.1,0.1)))

#model reporting
inlaRANEF = round(inla1$summary.random$id[,c(2,4,6)],2)
plot(1:103,inlaRANEF[,1], pch=1, ylim=c(min(inlaRANEF[,2]),max(inlaRANEF[,3])),main="Posterior E[b_i|Y] and 95% Credible Intervals", xlab="Cluster Index i", ylab="E[b_i|Y]")
points((1:103)+0.3,ranef(lmer1)$id[,1],pch=17)
for(i in 1:103)
{
  segments(x0=i, y0=inlaRANEF[i,2], y1=inlaRANEF[i,3])
}
legend("topleft", legend=c("INLA", "LMER"), pch=c(1,17))

# fixed effects
xtable(round(rbind(inla1$summary.fixed[,c(4,3,5)],1/sqrt(inla1$summary.hyperpar[,c(4,3,5)])),2))

#diagnostics
ordered = order(inla1$summary.random$id[,2]) #just order the posterior means of the b_i
x = qnorm(ppoints(103)) #need these
y = inla1$summary.random$id[ordered,2]
lower = inla1$summary.random$id[ordered,4]
upper = inla1$summary.random$id[ordered,6]
plot(x,y,ylim=c(min(lower),max(upper)), xlab="qnorm(ppoints(103))", ylab="Ordered E[b_i|Y]")
for(i in 1:103)
{
  segments(x0=x[i], y0=lower[i], y1=upper[i])
}
plot(inla1$summary.random$id[,3]^(-2), inla1$summary.random$id[,2],
     xlab="Precision = 1/Posterior Variance", ylab="Posterior Mean") #WHAT DOES THIS SHOW

#now for LOO
set.seed(1)
beta0s = matrix(0,nrow=103,ncol=3)
beta1s = beta0s
beta2s = beta0s
beta3s = beta0s
beta0sLMER = beta0s
beta1sLMER = beta0s
beta2sLMER = beta0s
beta3sLMER = beta0s
for(i in 1:103)
{
  tempDF = subset(na.omit(prob2DF),id!=i)
  inlaLOO <- inla(formula1, family="gaussian", data=tempDF, control.fixed=fixedprior, control.family = list(prior="loggamma",param=c(0.1,0.1)))
  lmerLOO <- lmer(y~trt+time+trt*time + (1 |id) ,data=tempDF)
  beta0s[i,] = as.numeric(inlaLOO$summary.fixed[1,c(1,3,5)])
  beta1s[i,] = as.numeric(inlaLOO$summary.fixed[2,c(1,3,5)])
  beta2s[i,] = as.numeric(inlaLOO$summary.fixed[3,c(1,3,5)])
  beta3s[i,] = as.numeric(inlaLOO$summary.fixed[4,c(1,3,5)])
  beta0sLMER[i,] = summary(lmerLOO)$coefficients[1,1]+c(0,-1,1)*1.96*summary(lmerLOO)$coefficients[1,2]
  beta1sLMER[i,] = summary(lmerLOO)$coefficients[2,1]+c(0,-1,1)*1.96*summary(lmerLOO)$coefficients[2,2]
  beta2sLMER[i,] = summary(lmerLOO)$coefficients[3,1]+c(0,-1,1)*1.96*summary(lmerLOO)$coefficients[3,2]
  beta3sLMER[i,] = summary(lmerLOO)$coefficients[4,1]+c(0,-1,1)*1.96*summary(lmerLOO)$coefficients[4,2]
}
plot(beta0s[,1],ylim=c(200,250),main="Intercept LOO Diagnostics")
points((1:103)+0.3,beta0sLMER[,1],pch=17)
for(i in 1:103)
{
  segments(x0=i, y0=beta0s[i,2], y1=beta0s[i,3])
  segments(x0=i+0.3, y0=beta0sLMER[i,2], y1=beta0sLMER[i,3])
}
legend("topleft", pch=c(1,17), legend=c("Bayes", "LMER"))

#beta1:
plot(beta1s[,1],ylim=c(-20,40),main="Beta1 LOO Diagnostics")
points((1:103)+0.3,beta1sLMER[,1],pch=17)
for(i in 1:103)
{
  segments(x0=i, y0=beta1s[i,2], y1=beta1s[i,3])
  segments(x0=i+0.3, y0=beta1sLMER[i,2], y1=beta1sLMER[i,3])
}
legend("topleft", pch=c(1,17), legend=c("Bayes", "LMER"))

#beta2:
plot(beta2s[,1],ylim=c(0,15),main="Beta2 LOO Diagnostics")
points((1:103)+0.5,beta2sLMER[,1],pch=17)
for(i in 1:103)
{
  segments(x0=i, y0=beta2s[i,2], y1=beta2s[i,3])
  segments(x0=i+0.5, y0=beta2sLMER[i,2], y1=beta2sLMER[i,3])
}
legend("topleft", pch=c(1,17), legend=c("Bayes", "LMER"))

#beta3
plot(beta3s[,1],ylim=c(-10,10),main="Beta3 LOO Diagnostics")
points((1:103)+0.3,beta3sLMER[,1],pch=17)
for(i in 1:103)
{
  segments(x0=i, y0=beta3s[i,2], y1=beta3s[i,3])
  segments(x0=i+0.3, y0=beta3sLMER[i,2], y1=beta3sLMER[i,3])
}
legend("topleft", pch=c(1,17), legend=c("Bayes", "LMER"))
###########
#problem 4#
###########

set.seed(572)
sb = 1
sy = 0.5
beta0 = 1
beta1 = -1
beta2 = 0.5
n = 10^3
repSize = 10^3
betas = matrix(0,nrow=repSize,ncol=3)
betasLMM = matrix(0,nrow=repSize, ncol=3)
count0 = 0
count1 = 0
count2 = 0
count0LMM = 0
count1LMM = 0
count2LMM = 0
for(k in 1:repSize)
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
      if(dat$Y[i+j] < -1)
      {
        missing$Y[((i+j+1):(i+2))] <- NA
        break
      }
    }
  }
  model3 = geem(Y~t+X,id=id,data=missing,corstr="exchangeable")
  #model2 = gee(Y~t+X, id=id, data=missing, na.action=na.omit, corstr="exchangeable")
  #model3 = geem(Y~t+X,id=id,data=missing,corstr="independence")
  modelLMM <- lmer(Y~t+X+(1|id), data=missing)
  betas[k,] = model3$beta
  betasLMM[k,] = fixef(modelLMM)
  len = 1.96*summary(model3)$se.robust
  lenLMM = 1.96*summary(modelLMM)$coefficients[,2]
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
  #####coverage intervals for LMM
  if(
    findInterval(beta0,fixef(modelLMM)[1] + c(-1,1)*lenLMM[1]) == 1
  )
  {
    count0LMM = count0LMM + 1
  }
  if(
    findInterval(beta1,fixef(modelLMM)[2] + c(-1,1)*lenLMM[2]) == 1
  )
  {
    count1LMM = count1LMM + 1
  }
  if(
    findInterval(beta2,fixef(modelLMM)[3] + c(-1,1)*lenLMM[3]) == 1
  )
  {
    count2LMM = count2LMM + 1
  }
}
hist(betas[,1],main=expression(paste("Histogram of ", hat(beta)[0], " (GEE)")),xlab=expression(hat(beta)[0]))
hist(betas[,2],main=expression(paste("Histogram of ", hat(beta)[1], " (GEE)")),xlab=expression(hat(beta)[1]))
hist(betas[,3],main=expression(paste("Histogram of ", hat(beta)[2], " (GEE)")),xlab=expression(hat(beta)[2]))
hist(betasLMM[,1],main=expression(paste("Histogram of ", hat(beta)[0], " (LMM)")),xlab=expression(hat(beta)[0]))
hist(betasLMM[,2],main=expression(paste("Histogram of ", hat(beta)[1], " (LMM)")),xlab=expression(hat(beta)[1]))
hist(betasLMM[,3],main=expression(paste("Histogram of ", hat(beta)[2], " (LMM)")),xlab=expression(hat(beta)[2]))
###########
#problem 5#
###########
ni = 3
n=10^3
p=1/4
theta0 = -2
theta1 = 1.5
beta1 = 0
beta2 = 1.5
sigb = 0.1

expit <- function(z)
{
  return(
    1/(1+exp(-z))
  )
}

reps=10^3
beta0 = -2.5
ps = matrix(NA, nrow=reps,ncol=2)
set.seed(571)
for(i in 1:reps)
{
  Z = rbinom(n, 1, prob=p)
  XList = lapply(expit(theta0+theta1*Z), function(p) rbinom(ni,2,p))
  b = rnorm(n,sd=sigb)
  Yprobs = expit(beta0+b+beta2*Z) #because beta1=0
  Y = lapply(Yprobs, function(p) rbinom(ni,1,p))
  prob5DF = data.frame(Y = unlist(Y), X = unlist(XList), Z=rep(Z,each=ni), id=rep(1:n, each=ni))
  tryCatch({
  m1 = lme(Y~X+Z, random = ~1|id, data=prob5DF)
  m2 = glmer(Y~X+Z + (1|id), family=binomial, data=prob5DF)
  ps[i,]=c(summary(m1)$tTable[2,5],summary(m2)$coeff[2,4])}, error = function(e){})
}
colMeans(na.omit(ps)<0.05)
hist(ps[,2],main="Histogram of p-values from GLMM, beta0=-2.5", xlab="p values")
#for an average cluster,
expit(beta0+beta2*p)


############
#JAGS CODE##
############
my.inits <- list( # set up initial values, for 2 chains;
   list( b0=rep(0,103), beta0=0, beta1=0, beta2=0, beta3=0, taub=1, tauY=1),
   list( b0=seq(0,1,l=103), beta0=-2, beta1=0.5, beta2=0, beta3=1, taub=0.5, tauY=1.5)
)
temp = na.omit(prob2DF)
prob2.long <- list( # data in list format
   y = temp$y,
  trt = temp$trt,
  time = temp$time,
  id = temp$id
)
jags = jags.model("p2bugscode.txt", data=prob2.long,
           inits=my.inits,
           n.chains=2)
mine = coda.samples(model = jags,
             variable.names = c("b0", "beta0","beta1","beta2", "beta3", "taub", "tauY"),
             n.iter = 200000)
samples.out.keep = window(mine, start=2001, end=200000) #burn first 2000
test = summary(mine)

indices = seq(150000,200000,by=500)

chain1 = mine[[1]]
fitted = matrix(0,nrow=447,ncol=length(indices))
resid = fitted
for(j in 1:length(indices))
{ #for each sample from the posterior,
  for(i in 1:447)
  { #for each element of the data frame, get the fitted and residuals
    fitted[i,j] = chain1[indices[j],temp[i,1]] + chain1[indices[j],104:107]%*%unlist(c(1,temp[i,c(4,2)],temp[i,4]*temp[i,2]))
    resid[i,j] = temp[i,3] - fitted[i]
  }
}
dim(fitted) = c(447*length(indices),1)
dim(resid) = c(447*length(indices),1)
plot(fitted,resid,col=rgb(0,0,0,alpha=0.1),main="Post. Resid. vs. Fitted for 100 Well-Spaced MCMC Samples")
lines(lowess(fitted,resid,iter=0),col="green")
abline(h=0)

plot(fitted(lmer1),residuals(lmer1,type="pearson"),ylim=c(-75,75),xlim=c(145,400),main="Residuals vs Fitted, LMM", ylab="Pearson Residuals (level 1)", xlab="Fitted Values")
lines(lowess(fitted(lmer1),residuals(lmer1,type="pearson")),col="green")
