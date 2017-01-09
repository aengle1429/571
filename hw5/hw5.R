library(gee)
library(geeM)
library(ggplot2)
library(reshape2)
library(xtable)
########################
#######Problem 1########
########################
data(Orthodont, package="nlme")
d4 <- Orthodont # using shorter names
d4$id <- d4$Subject
d4$male <- d4$Sex=="Male"
glm1 <- glm(distance~I(age-8)*male, data=d4)
glm2 <- glm(distance~age, data=d4)
gee1 <- geem(distance~age, id=id, data=d4, corstr="independence")
X = model.matrix(glm2)
A = sum((X[,2]-mean(X[,2]))^2)
B = 0
for(i in seq(1,length(d4$distance),by=4))
{
  B = B + sum((X[i:(i+3),2]-11)*(d4$distance[i:(i+3)]-X[i:(i+3),]%*%gee1$beta))^2
}

B/A^2

bhat <- sum((X[,2]-mean(X[,2]))*d4$distance)/A
#agrees with gee1

########################
#######Problem 2########
########################
set.seed(200)
n = 20 #20 chirrens
beta0 = 60
beta1 = 40 #english is harder
beta2 = 2 #age, Z
beta3 = 10 #math, X
beta4 = -10 #english, U
beta5 = 3
beta6 = -1
i = 1
data = data.frame()
repeat
{
  english = sample(1:n,size=(n/2),replace=F)
  math = sample(1:n,size=(n/2),replace=F)
  Z = sample(-5:5,size=n,replace=T)
  y1 = rep(beta0,n)
  y2 = rep(beta1,n)
  y1[math] = y1[math] + beta3 + Z[math]*beta5
  y2[english] = y2[english] + beta4 + Z[english]*beta6
  W = rep(0,n) #parental help
  help = sample(1:n, size=(n/2),replace=F)
  W[help] = 5 #no help
  y1 = y1 + beta2*Z +  rnorm(n,0,1) + W
  y2 = y2 + beta2*Z + rnorm(n,0,1) + W
  mathind = rep(0,n)
  mathind[math] = 1
  engind = rep(0,n)
  engind[english] = 1
  
  tempdata = data.frame(cbind(y1,y2)) #write y2 as percentage points
  start = 20*(i-1) + 1
  tempdata$id = factor(start:(start+19))
  sort = melt(tempdata,id.var="id")
  sort = sort[order(sort$id),]
  sort$Z = rep(Z,each=2)
  sort$X = rep(mathind,each=2)
  sort$U = rep(engind,each=2)
  data = rbind(data,sort)
  if(i >= 1000)
  {
    break
  }
  i = i + 1
}

data$U[seq(1,dim(data)[1],by=2)]=0 #indicator turn odd off
data$X[seq(2,dim(data)[1],by=2)]=0
data$int = rep(1,dim(data)[1])
data$int[seq(1,dim(data)[1],by=2)]=0 #fit separate intercepts
gee(value~int+Z+X+U+Z*X+Z*U,id=id, data=data, corstr="exchangeable")

########################
#######Problem 3########
########################
PI <- function(alpha,c,PE)
{
  return(
    1 - pnorm(
      (1-sqrt(c))*qnorm(1-alpha) + qnorm(1-PE)
              )
    )
}


q3 = matrix(0,nrow=105,ncol=2)
q3[,1] = rep(seq(0,1,0.05),5)
cs = c(0.1,0.3,0.5,0.7,0.9)
inds = (0:4)*21 + 1
j = 1
for(ind in inds)
{
  q3[ind:(ind+20),2] = PI(0.05,cs[j],seq(0,1,0.05))
  j = j + 1
}

df = data.frame(q3)
names(df) = c("PE", "PI")
df$ARE <- rep(as.factor(cs),each=21)
ggplot(df, aes(x=PE, y=PI,col=ARE)) + geom_point() + geom_line() + ggtitle(expression(paste("Power of Inefficient Working Covariance vs. Efficient Working Covariance at ", alpha, " = 0.05"))) +
  xlab("Power for Efficient Working Covariance") + ylab("Power for Inefficient Working Covariance")

#now do it for a new alpha
j = 1
for(ind in inds)
{
  q3[ind:(ind+20),2] = PI(1e-6,cs[j],seq(0,1,0.05))
  j = j + 1
}

df = data.frame(q3)
names(df) = c("PE", "PI")
df$ARE <- rep(as.factor(cs),each=21)
ggplot(df, aes(x=PE, y=PI,col=ARE)) + geom_point() + geom_line() + ggtitle(expression(paste("Power of Inefficient Working Covariance vs. Efficient Working Covariance at ", alpha, " = 1e-6"))) +
  xlab("Power for Efficient Working Covariance") + ylab("Power for Inefficient Working Covariance")


########################
#######Problem 4########
########################
expit <- function(x)
{
  return(
    exp(x)/(1+exp(x))
  )
}
#part (a)
set.seed(100)
ARE = matrix(0,nrow=6,ncol=3) #one row per sigma, one column per coefficient
m = 10^2
n = 20*m #clusters are integral multiples of 20.
ni = 6 #same for each clsuter
#clusters have same X covariate
x = (1:6-3)/3
#sigmas = seq(0.5,3,0.5)
sigmas = 3
beta0 = -2.5
beta1 = 1
beta2 = 1
beta = c(beta0,beta1,beta2)

X = rep(x,n)
for(j in 1:length(sigmas))
{
  b = sigmas[j]*rep(qnorm(ppoints(20)),each=(ni*m))
  #fixing b this way puts clusters with equal values of b_i next to one another

  ones = sample(1:n,size=(n/2),replace=F)
  z = c()
  for(i in 1:n)
  {
    if(i %in% ones)
    {
      z = append(z,rep(1,6))
    }
    else{
      z = append(z,rep(0,6))
    }
  }
  
  probs = expit(cbind(rep(1,(n*ni)), X, z)%*%beta + b)
  Y = rep(0,(n*ni))
  for(i in 1:(n*ni))
  {
    Y[i] = rbinom(1,1,probs[i])
  }
  id = rep(1:n,each=6) #cluster id
  dat = data.frame(cbind(Y,X,z,id,b)) #make the data frame
  dat2 = dat
  #COMMENT OUT FOR part(a)
  missing = sample(1:n, size=(n/2), replace=F) #which clusters to omit
  where = 6*(missing-1) + 1
  dat2 = dat[-c(where + 3 , where + 4, where + 5),]
  #Comment above... dont reseed
  geeIND = gee(Y~X+z,data=dat2, id=id, family=binomial, corstr="independence")
  geeEXCH = gee(Y~X+z,data=dat2,id=id, family=binomial, corstr="exchangeable")
  #ARE[j,]=summary(geeEXCH)$se.robust/summary(geeIND)$se.robust
  ARE[j,]=summary(geeEXCH)$coefficients[,4]^2/summary(geeIND)$coefficients[,4]^2
}

xtable(ARE)



