library(geeM)
library(arm)
library(ggplot2)
library(texreg)
library("lme4") # contains the (balanced) sleepstudy data
sleepstudy$DaysC <- with(sleepstudy, Days - mean(Days)) # centering
library("nlme")
lm0c <- lm( Reaction ~ DaysC, sleepstudy) # ignore clustering
lme1 <- lme(Reaction ~ DaysC, random=~1|Subject, data=sleepstudy)
lme2a <- lme(Reaction ~ DaysC, random=reStruct(~DaysC|Subject, pdClass="pdDiag"),
             data=sleepstudy) # diagonal G
lme2b <- lme(Reaction ~ DaysC, random=~DaysC|Subject,
             data=sleepstudy) # unrestricted G
yMarg <- fitted(lme2a, level=0)[1:10] #marginal means, just need to grab the first 10
#need Sigma_i = ziGzi^T + phi R_i now...
VarCorr(lme2a)
z = model.matrix(lme2a,data=sleepstudy)[1:10,] #need these rows for the z_i
G = diag(
  VarCorr(lme2a)[1:2,1]
)
#correlation DEFAULTS to null, so no correlation, so R_i = I
Sig = z%*%G%*%t(z) + diag(VarCorr(lme2a)[3,1],nrow=10,ncol=10)
ses = sqrt(diag(Sig))


upper = yMarg + 1.96*ses
lower = yMarg - 1.96*ses
plot(yMarg, ylim=c(100,500), main="Marignal Mean and 95% Confidence Intervals, Centered X")
for(i in 1:10)
{
  segments(x0=i,y0=lower[i],x1=i,y1=upper[i])
}

######now do it without centering, U is for uncentered!
lme2c <- lme(Reaction ~ Days, random=reStruct(~Days|Subject, pdClass="pdDiag"),
             data=sleepstudy) # diagonal G
yMargU <- fitted(lme2c, level=0)[1:10] #marginal means, just need to grab the first 10

zU = model.matrix(lme2c,data=sleepstudy)[1:10,] #need these rows for the z_i
GU = diag(
  VarCorr(lme2c)[1:2,1]
)
#correlation DEFAULTS to null, so no correlation, so R_i = I
SigU = zU%*%GU%*%t(zU) + diag(VarCorr(lme2c)[3,1],nrow=10,ncol=10)
sesU = sqrt(diag(SigU))


upperU = yMargU + 1.96*sesU
lowerU = yMargU - 1.96*sesU

plot(yMargU, ylim=c(100,500), main="Marignal Mean and 95% Confidence Intervals, Uncentered X")
for(i in 1:10)
{
  segments(x0=i,y0=lowerU[i],x1=i,y1=upperU[i])
}  
###########
#problem 2#
###########
n = 10^3
ni = 5
sigb = 1
sigy = 6
b = rnorm(1000,mean=0,sd=sigb)
do.one <- function(n)
{
  y = c()
  for(i in 1:n)
  {
    yi <- rnorm(ni, mean=b[i], sd = sigy) 
    y = append(y,yi)
  }
  mydat = data.frame(y=y,id=rep(1:n,each=ni))
  m1 <- lmer(y ~ 1|id, data=mydat)
  si = se.ranef(m1)$id #same for all, just grab first
  bi = ranef(m1)$id
  lower = bi - 1.96*si
  upper = bi + 1.96*si
  return(mean((lower < b) & (b < upper)))
}

P = replicate(10^4, do.one(10^3))
##################
#problems 3 and 4#
##################
setwd("~/Dropbox/UW2015-2016/Win2016/571/hw8")
setwd("C:\\Users\\aengl_000\\Dropbox\\UW2015-2016\\Win2016\\571\\hw8")
dat <- read.table("xerop.csv", sep=",", header=T)

dat$cos <- cos(0.5*pi*(dat$time+1))
dat$sin <- sin(0.5*pi*(dat$time+1))

m1 <- geem(respinfect ~ sex + ht.for.age  + cos + sin + xerop + poly(age,degree=2,raw=T), family=binomial, corstr="ar1", data=dat, id=id)

m2 <- glmer(respinfect ~ (1|id) + sex + ht.for.age  + cos + sin + xerop +  poly(age,degree=2,raw=T), family=binomial, data=dat, start=list(fixef = m1$beta))
m3 <- glmer(respinfect ~ (1|id) + sex + ht.for.age  + cos + sin  +  poly(age,degree=2,raw=T), family=binomial, data=dat)
xtable(anova(m3,m2))
############
#LOO diag###
############
# ans = c()
# for(ids in unique(dat$id))
# {
#    leftOut = subset(dat,id!=ids)
#    tempMod <- glmer(respinfect ~ (1|id) + sex + ht.for.age  + cos + sin + xerop + poly(age, degree=2,raw=T), family=binomial, data=leftOut, start=list(fixef = fixef(m2)))
#    se = sqrt(diag(vcov(tempMod)))
#    ans = rbind(ans,cbind(Est = fixef(tempMod), LL = fixef(tempMod) - 1.96 * se, UL = fixef(tempMod) + 1.96 * se, id=0:7))
#    cat("id is ", ids, "\n")
# }
# LOOdata = data.frame(ans, cluster = rep(1:275, each=7))
# 
# 
# plotbeta0 = subset(data, id==0)
# ggplot(plotbeta0, aes(x = cluster, y = Est)) +
#   geom_errorbar(aes(ymin=LL, ymax=UL), width=.1) + 
#   geom_point(size=1)+ggtitle("beta0: intercept")
# 
# plotbeta1 = subset(data, id==1)
# ggplot(plotbeta1, aes(x = cluster, y = Est)) +
#   geom_errorbar(aes(ymin=LL, ymax=UL), width=.1) + 
#   geom_point(size=1)+ggtitle("beta1")
# 
# plotbeta2 = subset(data, id==2)
# ggplot(plotbeta2, aes(x = cluster, y = Est)) +
#   geom_errorbar(aes(ymin=LL, ymax=UL), width=.1) + 
#   geom_point(size=1)+ggtitle("beta2")
# 
# plotbeta3 = subset(data, id==3)
# ggplot(plotbeta3, aes(x = cluster, y = Est)) +
#   geom_errorbar(aes(ymin=LL, ymax=UL), width=.1) + 
#   geom_point(size=1)+ggtitle("beta3")
# 
# plotbeta4 = subset(data, id==4)
# ggplot(plotbeta4, aes(x = cluster, y = Est)) +
#   geom_errorbar(aes(ymin=LL, ymax=UL), width=.1) + 
#   geom_point(size=1)+ggtitle("beta4")
# 
# plotbeta5 = subset(data, id==5)
# ggplot(plotbeta5, aes(x = cluster, y = Est)) +
#   geom_errorbar(aes(ymin=LL, ymax=UL), width=.1) + 
#   geom_point(size=1)+ggtitle("beta5")
# 
# plotbeta6 = subset(data, id==6)
# ggplot(plotbeta6, aes(x = cluster, y = Est)) +
#   geom_errorbar(aes(ymin=LL, ymax=UL), width=.1) + 
#   geom_point(size=1)+ggtitle("beta6")


#subject-level fitted values against subject-level pearson residuals (level 1 in both cases)
plot(fitted(m2), residuals(m2, type="pearson"),ylab="PearsonResiduals",xlab="Fitted Values", main="Pearson Residuals vs. Fitted Values")
lines(lowess(fitted(m2),residuals(m2,type="pearson"),iter=0),col="red")
abline(h=0, lty=2)

#Pearson Residuals vs. sex, ht.for.age, cos, sin, xerop, age
plot(dat$sex,residuals(m2,type="pearson"),ylab="PearsonResiduals",xlab="sex", main="Pearson Residuals vs. Sex")
lines(lowess(dat$sex,residuals(m2,type="pearson"),iter=0),col="blue",lwd=2)
abline(h=0, lty=2)

plot(dat$ht.for.age,residuals(m2,type="pearson"),ylab="PearsonResiduals",xlab="ht.for.age", main="Pearson Residuals vs. Ht.For.Age")
lines(lowess(dat$ht.for.age,residuals(m2,type="pearson"),iter=0),col="blue",lwd=2)
abline(h=0, lty=2)

plot(dat$cos,residuals(m2,type="pearson"),ylab="PearsonResiduals",xlab="cos", main="Pearson Residuals vs. Cosine")
lines(lowess(dat$cos,residuals(m2,type="pearson"),iter=0),col="blue",lwd=2)
abline(h=0,lty=2)

plot(dat$sin,residuals(m2,type="pearson"),ylab="PearsonResiduals",xlab="sin",main="Pearson Residuals vs. Sine")
lines(lowess(dat$sin,residuals(m2,type="pearson"),iter=0),col="blue",lwd=2)
abline(h=0, lty=2)

plot(dat$xerop,residuals(m2,type="pearson"),ylab="PearsonResiduals",xlab="xerop", main="Pearson Residuals vs. Dry Eye Presence")
lines(lowess(dat$xerop,residuals(m2,type="pearson"),iter=0),col="blue",lwd=2)
abline(h=0, lty=2)

plot(dat$age,residuals(m2,type="pearson"),ylab="PearsonResiduals",xlab="age", main="Pearson Residuals vs. Age")
lines(lowess(dat$age,residuals(m2,type="pearson"),iter=0),col="blue",lwd=2)
abline(h=0, lty=2)

plot(dat$age^2,residuals(m2,type="pearson"),ylab="PearsonResiduals",xlab="age", main="Pearson Residuals vs. Age^2")
lines(lowess(dat$age^2,residuals(m2,type="pearson"),iter=0),col="blue",lwd=2)
abline(h=0, lty=2)


###variance diagnostic
# plot(fitted(m2), residuals(m2, type="pearson")^2,log="y", main="Sq Pearson Residuals vs. Fitted Values")
# lines(lowess(fitted(m2),residuals(m2,type="pearson")^2,iter=0),col="red")
# abline(h=1, lty=2)
# 
# #Pearson Residuals vs. sex, ht.for.age, cos, sin, xerop, age
# plot(dat$sex,residuals(m2,type="pearson")^2,ylab="PearsonResiduals",xlab="sex",log="y", main="Sq Pearson Residuals vs. Sex")
# lines(lowess(dat$sex,residuals(m2,type="pearson")^2,iter=0),col="blue",lwd=2)
# abline(h=1)
# 
# plot(dat$ht.for.age,residuals(m2,type="pearson")^2,ylab="PearsonResiduals",xlab="ht.for.age", log="y", main="Sq Pearson Residuals vs. Ht.For.Age")
# lines(lowess(dat$ht.for.age,residuals(m2,type="pearson")^2,iter=0),col="blue",lwd=2)
# abline(h=1)
# 
# plot(dat$cos,residuals(m2,type="pearson")^2,ylab="PearsonResiduals",xlab="cos", log="y", main="Sq Pearson Residuals vs. Cosine")
# lines(lowess(dat$cos,residuals(m2,type="pearson")^2,iter=0),col="blue",lwd=2)
# abline(h=1)
# 
# plot(dat$sin,residuals(m2,type="pearson")^2,ylab="PearsonResiduals",xlab="sin", log="y", main="Sq Pearson Residuals vs. Sine")
# lines(lowess(dat$sin,residuals(m2,type="pearson")^2,iter=0),col="blue",lwd=2)
# abline(h=1)
# 
# plot(dat$xerop,residuals(m2,type="pearson")^2,ylab="PearsonResiduals",xlab="xerop", log="y", main="Sq Pearson Residuals vs. Dry Eye")
# lines(lowess(dat$xerop,residuals(m2,type="pearson")^2,iter=0),col="blue",lwd=2)
# abline(h=1)
# 
# plot(dat$age,residuals(m2,type="pearson")^2,ylab="PearsonResiduals",xlab="age", log="y", main="Sq Pearson Residuals vs Age")
# lines(lowess(dat$age,residuals(m2,type="pearson")^2,iter=0),col="blue",lwd=2)
# abline(h=1)
# 
# #qqplot
# qqnorm((ranef(m2)$id[,1]))
# qqline((ranef(m2)$id[,1]))
# hist(ranef(m2)$id[,1])
