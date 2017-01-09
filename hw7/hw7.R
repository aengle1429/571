setwd("~/Dropbox/UW2015-2016/Win2016/571/hw7")

dat <- read.table("dyestuff.txt", header = T)

wide = reshape(dat, timevar = "rep", idvar = "batch", direction = "wide")
wide <- wide[,-1]

m = 5
n = 6
mu = mean(dat[, 1]) #grand mean
Yibar <- apply(wide, 1, mean) #cluster means

f = m/(n-1)*sum((Yibar - mu)^2)*n*(m-1)/sum((wide-Yibar)^2)
1-pf(f,n-1, n*(m - 1))


library(nlme)
m1 <- lm(dist~speed+I(speed^2), data=cars)
dat = m1$residuals
dat2 = data.frame(dat,rep(c(1,0),each=25),1:50)
names(dat2) = c('y', 'z', 'id')
t.test(formula=y~z, data=dat2)
summary(fm1 <- lme(fixed=y ~ z, random=~z|id, data=dat2,method='ML'))
