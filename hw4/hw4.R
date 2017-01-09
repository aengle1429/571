library(nlme)
library(ggplot2)
library(geeM)
library(reshape2)
data(Orthodont, package="nlme")
d4 <- Orthodont # using shorter names
d4$id <- d4$Subject
d4$male <- d4$Sex=="Male"
glm1 <- glm(distance~I(age-8)*male, data=d4)
X = model.matrix(glm1)


coplot(distance~age|Subject,data=d4,show.given=F,
       col=c(rep("blue",64),rep("magenta",44)),
       xlim=c(8,14),ylim=c(15,32),
       panel=function(x,y,col,...){points(x,y,col=col)
       lines(x,y,col=col)})
# coplot(distance~age|Subject,data=d4,show.given=F,
#        col=c(rep("blue",48),rep("black",4),rep("blue",12),rep("magenta",8),rep("green",4),rep("magenta",32)),
#        xlim=c(8,14),ylim=c(15,32),
#        panel=function(x,y,col,...){points(x,y,col=col)
#          lines(x,y,col=col)})

dist = matrix(d4$distance,nrow=4,ncol=27,byrow=FALSE)
matplot(d4$age[1:4],dist,type="l",col=c(rep("blue",16),rep("magenta",11)),lty=1,
        xlim=c(8,14),xlab="Age", ylab="Distance (mm)", main="Spaghetti Plot Colored by Gender")
#MALE SUBJECT 13
#his start index is
#(13-1)*4+1
#[1] 49
#end is 53

#19 subject, female #3
#73-76

###first bootstrap
B=1000
fbetaBoot <- matrix(0,nrow=B,ncol=4)

d4blind <- d4
d4blind['blind']=rep(seq(1:27),each=4) #just give me the 27 kids
set.seed(1)
wide <- reshape(d4blind,timevar="age",
                idvar=c("Subject","Sex","id","male","blind"),dir="wide")
for(b in 1:B)
{
  vars = sample(1:27,size=27,replace=TRUE) #indices to use in bootstrap GEE
  boot = wide[vars,] #grab those kids from the WIDE data frame
  long = melt(boot,id.vars=c("Subject","Sex","id","male","blind"),
              measure.vars=c("distance.8","distance.10","distance.12","distance.14"))  
  long = long[as.vector(matrix(1:108,nrow=4,byrow=T)),]
  long['age'] = rep(c(8,10,12,14),27)
  long['bootID'] = rep(1:27,each=4) #have to include this otherwise multiple kids appearing won't be distinguished in GEE
  gee1temp <- geem(value~I(age-8)*male, id=bootID, data=long, corstr="ar1")
  fbetaBoot[b,] = gee1temp$beta
}

bootDF <- data.frame(as.vector(fbetaBoot),as.factor(rep(1:4,each=1000)))
names(bootDF) <- c("beta","index")

ggplot(bootDF, aes(x=beta,fill=index)) +
  geom_histogram(data=subset(bootDF, index=="1"), binwidth=0.25) +
  geom_histogram(data=subset(bootDF, index=="2"), binwidth=.025) +
  geom_histogram(data=subset(bootDF, index=="3"), binwidth=0.25) +
  geom_histogram(data=subset(bootDF, index=="4"), binwidth=.05) +
  facet_wrap("index", scales="free")


#now do it for Exchangeable...
fbetaBootE<- matrix(0,nrow=B,ncol=4)
for(b in 1:B)
{
  vars = sample(1:27,size=27,replace=TRUE) #indices to use in bootstrap GEE
  boot = wide[vars,] #grab those kids from the WIDE data frame
  long = melt(boot,id.vars=c("Subject","Sex","id","male","blind"),
              measure.vars=c("distance.8","distance.10","distance.12","distance.14"))  
  long = long[as.vector(matrix(1:108,nrow=4,byrow=T)),]
  long['age'] = rep(c(8,10,12,14),27)
  long['bootID'] = rep(1:27,each=4) #have to include this otherwise multiple kids appearing won't be distinguished in GEE
  gee1temp <- geem(value~I(age-8)*male, id=bootID, data=long, corstr="exchangeable")
  fbetaBootE[b,] = gee1temp$beta
}

bootDF2 <- data.frame(as.vector(fbetaBootE),as.factor(rep(1:4,each=1000)))
names(bootDF2) <- c("beta","index")

ggplot(bootDF2, aes(x=beta,fill=index)) +
  geom_histogram(data=subset(bootDF2, index=="1"), binwidth=0.25) +
  geom_histogram(data=subset(bootDF2, index=="2"), binwidth=.025) +
  geom_histogram(data=subset(bootDF2, index=="3"), binwidth=0.25) +
  geom_histogram(data=subset(bootDF2, index=="4"), binwidth=.05) +
  facet_wrap("index", scales="free")

#first ar1, diag for diagnostic
gee1 <- geem(distance~I(age-8)*male, id=id, data=d4, corstr="ar1")
phiAR1 <- gee1$phi
betaAR1 <- gee1$beta
d4DiagAR1 = d4
d4DiagAR1['cluster13'] = c(rep("blue",64),rep("magenta",44))
d4DiagAR1['resid'] = phiAR1^(-1/2)*(d4$distance - X%*%betaAR1)
d4DiagAR1['fitted'] = X%*%betaAR1
plot(d4DiagAR1$fitted,d4DiagAR1$resid,xlab="Fitted Values (Distance)",
     ylab="Pearson Residuals", 
     main="Pearson Residuals vs Fitted Values With Lowess iter=0",
     col=d4DiagAR1$cluster13,pch=c(rep(20,48),rep(3,4),rep(20,12),rep(20,8),rep(5,4),rep(20,32)))
lines(lowess(fitted.values(gee1),d4DiagAR1$resid,iter=0),col="green")

plot(rep(1:27,each=4),d4DiagAR1$resid,xlab="Child",
     ylab="Pearson Residuals", 
     main="Pearson Residuals vs Cluster Index With Lowess iter=0",
     col=d4DiagAR1$cluster13,pch=c(rep(20,48),rep(3,4),rep(20,12),rep(20,8),rep(5,4),rep(20,32)))
lines(lowess(rep(1:27,each=4),d4DiagAR1$resid,iter=0),col="green")


# unitLagAR1 = cbind(d4DiagAR1$resid[-108],d4DiagAR1$resid[-1])
# unitLagAR1 = unitLagAR1[-seq(4,104,4),]
# plot(unitLagAR1)
# abline(a=0,b=gee1$alpha)

plot(d4DiagAR1$fitted,d4DiagAR1$resid^2,xlab="Fitted Values (Distance)",
     ylab="Squared Pearson Residuals", 
     main="Squared P. Residuals vs Fitted Values With Lowess iter=0",
     col=d4DiagAR1$cluster13,
     pch=c(rep(20,48),rep(3,4),rep(20,12),rep(20,8),rep(5,4),rep(20,32)),
     log="y")
lines(lowess(d4DiagAR1$fitted,d4DiagAR1$resid^2,iter=0),col="green")
abline(h=1)

#semivariogram
semialpha1 = matrix(0,nrow=81,ncol=2) 
semialpha2 = matrix(0,nrow=81,ncol=2)
semialpha3 = matrix(0,nrow=81,ncol=2)
count = 1
for(i in seq(1,105,by=4))
{
  semialpha1[count,] = c(d4DiagAR1$resid[i],d4DiagAR1$resid[i+1])
  semialpha1[count+1,] = c(d4DiagAR1$resid[i+1],d4DiagAR1$resid[i+2])
  semialpha1[count+2,] = c(d4DiagAR1$resid[i+2],d4DiagAR1$resid[i+3])
  semialpha2[count,] = c(d4DiagAR1$resid[i],d4DiagAR1$resid[i+2])
  semialpha2[count+1,] = c(d4DiagAR1$resid[i+1],d4DiagAR1$resid[i+3])
  semialpha3[count,] = c(d4DiagAR1$resid[i],d4DiagAR1$resid[i+3])
  count = count + 3
}
semialpha2 = semialpha2[which(rowSums(semialpha2)!=0),]
semialpha3 = semialpha3[which(rowSums(semialpha3)!=0),]


semialpha1Est = 0.5*mean((semialpha1[,1]-semialpha1[,2])^2)
semialpha2Est = 0.5*mean((semialpha2[,1]-semialpha2[,2])^2)
semialpha3Est = 0.5*mean((semialpha3[,1]-semialpha3[,2])^2)

fitted = 1-c(gee1$alpha,gee1$alpha^2,gee1$alpha^3)
plot(fitted,c(semialpha1Est,semialpha2Est,semialpha3Est),xlim=c(0,1),ylim=c(0,1.2),
     main="Semivariogram for Working AR(1) Correlation", xlab="1-Fitted Element of Ri", ylab="0.5*mean(Delta(P.Resid)^2)")
abline(0,1)

#exchangeable
gee2 <- geem(distance~I(age-8)*male, id=id, data=d4, corstr="exchangeable")
phiExch <- gee2$phi
betaExch <- gee2$beta
d4DiagExch = d4
d4DiagExch['cluster13'] = c(rep("blue",64),rep("magenta",44))
d4DiagExch['resid'] = phiExch^(-1/2)*(d4$distance - X%*%betaExch)
d4DiagExch['fitted'] = X%*%betaAR1
plot(d4DiagExch$fitted,d4DiagExch$resid,xlab="Fitted Values (Distance)",
     ylab="Pearson Residuals", 
     main="Pearson Residuals vs Fitted Values With Lowess iter=0",
     col=d4DiagExch$cluster13,pch=c(rep(20,48),rep(3,4),rep(20,12),rep(20,8),rep(5,4),rep(20,32)))
lines(lowess(d4DiagExch$fitted,d4DiagExch$resid,iter=0),col="green")

plot(rep(1:27,each=4),d4DiagExch$resid,xlab="Child",
     ylab="Pearson Residuals", 
     main="Pearson Residuals vs Cluster Index With Lowess iter=0",
     col=d4DiagExch$cluster13,pch=c(rep(20,48),rep(3,4),rep(20,12),rep(20,8),rep(5,4),rep(20,32)))
lines(lowess(rep(1:27,each=4),d4DiagExch$resid,iter=0),col="green")


unitLagExch = cbind(d4DiagExch$resid[-108],d4DiagExch$resid[-1])
unitLagExch = unitLagExch[-seq(4,104,4),]
0.5*mean((unitLagExch[,1]-unitLagExch[,2])^2)
plot(unitLagExch)
abline(a=0,b=gee2$alpha)




#variance
plot(d4DiagExch$fitted,d4DiagExch$resid^2,xlab="Fitted Values (Distance)",
     ylab="Squared Pearson Residuals", 
     main="Squared P. Residuals vs Fitted Values With Lowess iter=0",
     log="y", col=c(rep("blue",64),rep("magenta",44)),pch=c(rep(20,48),rep(3,4),rep(20,12),rep(20,8),rep(5,4),rep(20,32)))
lines(lowess(d4DiagExch$fitted,d4DiagExch$resid^2,iter=0),col="green")
abline(h=1)

#LOO bootstraps
B=1000
beta0 <- matrix(0,nrow=B,ncol=27)
beta1 <- matrix(0,nrow=B,ncol=27)
beta2 <- matrix(0,nrow=B,ncol=27)
beta3 <- matrix(0,nrow=B,ncol=27)
d4blind <- d4
d4blind['blind']=rep(seq(1:27),each=4) #just give me the 27 kids
for(i in 1:27)
{
  start = (i-1)*4+1 #grab the starting index to leave one out
  d4temp <- d4blind[-c(start:(start+3)),] #give me the data frame w/o that kiddo
  wide <- reshape(d4temp,timevar="age",
                  idvar=c("Subject","Sex","id","male","blind"),dir="wide")
  for(b in 1:B)
  {
    vars = sample(1:26,size=26,replace=TRUE) #indices to use in bootstrap GEE
    boot = wide[vars,] #grab those kids from the WIDE data frame
    long = melt(boot,id.vars=c("Subject","Sex","id","male","blind"),
         measure.vars=c("distance.8","distance.10","distance.12","distance.14"))
    long = long[as.vector(matrix(1:104,nrow=4,byrow=T)),]
    long['age'] = rep(c(8,10,12,14),26)
    long['bootID'] = rep(1:26,each=4) #have to include this otherwise multiple kids appearing won't be distinguished in GEE
    gee1temp <- geem(value~I(age-8)*male, id=bootID, data=long, corstr="ar1")
    beta0[b,i] = gee1temp$beta[1]
    beta1[b,i] = gee1temp$beta[2]
    beta2[b,i] = gee1temp$beta[3]
    beta3[b,i] = gee1temp$beta[4]
  }
}


#plot bootstrap histograms for each parameter
beta1LOOdf <- data.frame(as.vector(beta1),as.factor(rep(1:27,each=1000)))
names(beta1LOOdf) <- c("beta1", "leftOut")
ggplot(beta1LOOdf, aes(x=leftOut, y=beta1)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data AR(1)") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[1])))

beta0LOOdf <- data.frame(as.vector(beta0),as.factor(rep(1:27,each=1000)))
names(beta0LOOdf) <- c("beta0", "leftOut")
ggplot(beta0LOOdf, aes(x=leftOut, y=beta0)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data AR(1)") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[0])))

beta2LOOdf <- data.frame(as.vector(beta2),as.factor(rep(1:27,each=1000)))
names(beta2LOOdf) <- c("beta2", "leftOut")
ggplot(beta2LOOdf, aes(x=leftOut, y=beta2)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data AR(1)") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[2])))

beta3LOOdf <- data.frame(as.vector(beta3),as.factor(rep(1:27,each=1000)))
names(beta3LOOdf) <- c("beta3", "leftOut")
ggplot(beta3LOOdf, aes(x=leftOut, y=beta3)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data AR(1)") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[3])))


#use GEE to get intervals with the cluster left out
d4blind <- d4
d4blind['blind']=rep(seq(1:27),each=4) #just give me the 27 kids
betas <- matrix(0,nrow=4,ncol=27)
CIbeta0s <- matrix(0,nrow=2,ncol=27)
CIbeta1s <- matrix(0,nrow=2,ncol=27)
CIbeta2s <- matrix(0,nrow=2,ncol=27)
CIbeta3s <- matrix(0,nrow=2,ncol=27)
for(i in 1:27)
{
  start = (i-1)*4+1 #grab the starting index to leave one out
  d4temp <- d4blind[-c(start:(start+3)),] #give me the data frame w/o that kiddo
  geetemp1 <- geem(distance~I(age-8)*male, id=id, data=d4temp, corstr="ar1")
  betas[,i] = summary(geetemp1)$beta
  CIbeta0s[,i] = summary(geetemp1)$se.robust[1]*1.96*c(-1,1)+betas[1,i]
  CIbeta1s[,i] = summary(geetemp1)$se.robust[2]*1.96*c(-1,1)+betas[2,i]
  CIbeta2s[,i] = summary(geetemp1)$se.robust[3]*1.96*c(-1,1)+betas[3,i]
  CIbeta3s[,i] = summary(geetemp1)$se.robust[4]*1.96*c(-1,1)+betas[4,i]
}
plot(1:27,betas[1,],ylim=c(19,23),col="red",pch=20,
     main="Intercept Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[0], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta0s[1,i],x1=i,y1=CIbeta0s[2,i])
}
plot(1:27,betas[2,],ylim=c(0.3,0.7),col="red",pch=20,
 main="Slope Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[1], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta1s[1,i],x1=i,y1=CIbeta1s[2,i])
}
plot(1:27,betas[3,],ylim=c(-0.5,3.5),col="red",pch=20,
     main="Slope Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[2], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta2s[1,i],x1=i,y1=CIbeta2s[2,i])
}
plot(1:27,betas[4,],col="red",ylim=c(0,0.6),pch=20,
     main="Slope Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[3], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta3s[1,i],x1=i,y1=CIbeta3s[2,i])
}


####now do all that boot LOO for exchangeable
B=1000
beta0 <- matrix(0,nrow=B,ncol=27)
beta1 <- matrix(0,nrow=B,ncol=27)
beta2 <- matrix(0,nrow=B,ncol=27)
beta3 <- matrix(0,nrow=B,ncol=27)
d4blind <- d4
d4blind['blind']=rep(seq(1:27),each=4) #just give me the 27 kids
set.seed(1)
for(i in 1:27)
{
  start = (i-1)*4+1 #grab the starting index to leave one out
  d4temp <- d4blind[-c(start:(start+3)),] #give me the data frame w/o that kiddo
  wide <- reshape(d4temp,timevar="age",
                  idvar=c("Subject","Sex","id","male","blind"),dir="wide")
  for(b in 1:B)
  {
    vars = sample(1:26,size=26,replace=TRUE) #indices to use in bootstrap GEE
    boot = wide[vars,] #grab those kids from the WIDE data frame
    long = melt(boot,id.vars=c("Subject","Sex","id","male","blind"),
                measure.vars=c("distance.8","distance.10","distance.12","distance.14"))
    long = long[as.vector(matrix(1:104,nrow=4,byrow=T)),]
    long['age'] = rep(c(8,10,12,14),26)
    long['bootID'] = rep(1:26,each=4) #have to include this otherwise multiple kids appearing won't be distinguished in GEE
    gee1temp <- geem(value~I(age-8)*male, id=bootID, data=long, corstr="exchangeable")
    beta0[b,i] = gee1temp$beta[1]
    beta1[b,i] = gee1temp$beta[2]
    beta2[b,i] = gee1temp$beta[3]
    beta3[b,i] = gee1temp$beta[4]
  }
}


#plot bootstrap histograms for each parameter
beta1LOOdf <- data.frame(as.vector(beta1),as.factor(rep(1:27,each=1000)))
names(beta1LOOdf) <- c("beta1", "leftOut")
ggplot(beta1LOOdf, aes(x=leftOut, y=beta1)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data Exchangeable") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[1])))

beta0LOOdf <- data.frame(as.vector(beta0),as.factor(rep(1:27,each=1000)))
names(beta0LOOdf) <- c("beta0", "leftOut")
ggplot(beta0LOOdf, aes(x=leftOut, y=beta0)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data Exchangeable") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[0])))

beta2LOOdf <- data.frame(as.vector(beta2),as.factor(rep(1:27,each=1000)))
names(beta2LOOdf) <- c("beta2", "leftOut")
ggplot(beta2LOOdf, aes(x=leftOut, y=beta2)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data Exchangeable") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[2])))

beta3LOOdf <- data.frame(as.vector(beta3),as.factor(rep(1:27,each=1000)))
names(beta3LOOdf) <- c("beta3", "leftOut")
ggplot(beta3LOOdf, aes(x=leftOut, y=beta3)) + geom_boxplot() + ggtitle("Deletion Diagnostics for Dental Data Exchangeable") + xlab("Child Not Appearing") + ylab(expression(paste("Values of ", hat(beta)[3])))


#use GEE to get intervals with the cluster left out
d4blind <- d4
d4blind['blind']=rep(seq(1:27),each=4) #just give me the 27 kids
betas <- matrix(0,nrow=4,ncol=27)
CIbeta0s <- matrix(0,nrow=2,ncol=27)
CIbeta1s <- matrix(0,nrow=2,ncol=27)
CIbeta2s <- matrix(0,nrow=2,ncol=27)
CIbeta3s <- matrix(0,nrow=2,ncol=27)
for(i in 1:27)
{
  start = (i-1)*4+1 #grab the starting index to leave one out
  d4temp <- d4blind[-c(start:(start+3)),] #give me the data frame w/o that kiddo
  geetemp1 <- geem(distance~I(age-8)*male, id=id, data=d4temp, corstr="exchangeable")
  betas[,i] = summary(geetemp1)$beta
  CIbeta0s[,i] = summary(geetemp1)$se.robust[1]*1.96*c(-1,1)+betas[1,i]
  CIbeta1s[,i] = summary(geetemp1)$se.robust[2]*1.96*c(-1,1)+betas[2,i]
  CIbeta2s[,i] = summary(geetemp1)$se.robust[3]*1.96*c(-1,1)+betas[3,i]
  CIbeta3s[,i] = summary(geetemp1)$se.robust[4]*1.96*c(-1,1)+betas[4,i]
}
plot(1:27,betas[1,],ylim=c(19,23),col="red",pch=20,
     main="Intercept Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[0], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta0s[1,i],x1=i,y1=CIbeta0s[2,i])
}
plot(1:27,betas[2,],ylim=c(0.3,0.7),col="red",pch=20,
     main="Slope Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[1], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta1s[1,i],x1=i,y1=CIbeta1s[2,i])
}
plot(1:27,betas[3,],ylim=c(-0.5,3.5),col="red",pch=20,
     main="Slope Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[2], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta2s[1,i],x1=i,y1=CIbeta2s[2,i])
}
plot(1:27,betas[4,],col="red",ylim=c(0,0.6),pch=20,
     main="Slope Estimate, Leaving out Cluster with 95% GEE Interval", xlab="Cluster Omitted", ylab=expression(paste(hat(beta)[3], "with Confidence Interval")))
for(i in 1:27)
{
  segments(x0=i,y0=CIbeta3s[1,i],x1=i,y1=CIbeta3s[2,i])
}