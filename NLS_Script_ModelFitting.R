#  -------------------------------------------------------------------------------
# Script and models used in model fitting procedure. 

# Based on figures originally presented in Carstensen and Weydmann 2012.

#  -------------------------------------------------------------------------------

library("nlme")

#pdf("models.pdf")
# Set plotting area
par(mfrow=c(2,3), mar=c(1,2,3,3), bty="n")


# --------------------------------------------------------------
# Null model (i.e mean) 
# --------------------------------------------------------------

x= seq(1,30, 1)
y = 10 + rnorm(length(x), 0, 5) 
y = (y-min(y))/(max(y)- min(y))
plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
axis(side=2, labels=F)
fit.gnls<-gnls(y~b0,start=c(b0=0))

lines(x, predict(fit.gnls))
lines(x, predict(fit.gnls)+sd(y), lty=2)
lines(x, predict(fit.gnls)-sd(y), lty=2)
text(2, 1.45, "A", cex=1.2 )
text(3, 1.45, "Constant mean", pos=4, cex=0.8)
text(3, 1.25, "Constant variance", pos=4, cex=0.8)



# --------------------------------------------------------------
# Constant mean with increasing variance.
# --------------------------------------------------------------

x= seq(1,30, 1)
y = 10 + (x/10)*rnorm(length(x), 0, 5) 
y = (y-min(y))/(max(y)- min(y))

plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
axis(side=2, labels=F)

fit.gnls<-gnls(y~b0,start=c(b0=0), weights =varFixed(~x))
sig = attr(residuals(fit.gnls), "std")
lines(x, predict(fit.gnls))
lines(x, predict(fit.gnls)+sig, lty=2)
lines(x, predict(fit.gnls)-sig, lty=2) 

text(2, 1.45, "B", cex=1.2 )
text(3, 1.45, "Constant mean", pos=4, cex=0.8)
text(3, 1.25, "Increasing std.error", pos=4, cex=0.8)





# --------------------------------------------------------------
# Linear model.
# --------------------------------------------------------------

x= seq(1,30, 1)
y = 10 + (2*x) + rnorm(length(x), 0, 5) 
y = (y-min(y))/(max(y)- min(y))

plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
axis(side=2, labels=F)
fit.gnls<-gnls(y~b0+b1*(x),start=c(b0=0,b1=0) )
sig = summary(fit.gnls)$sigma # standard deviation of the residuals (called standard error in R)
lines(x, predict(fit.gnls))
lines(x, predict(fit.gnls)+sig, lty=2)
lines(x, predict(fit.gnls)-sig, lty=2) 

text(2, 1.45, "C", cex=1.2 )
text(3, 1.45, "Linear relationship", pos=4, cex=0.8)
text(3, 1.25, "Constant std.error", pos=4, cex=0.8)



# --------------------------------------------------------------
# Segmented linear model.
# --------------------------------------------------------------

x= seq(1,30, 1)
k= 15; d=4
y = ifelse(x-k<= 0, x/10, x/10 + d*(x-k)) + rnorm(length(x), 0, 5) 
y = (y-min(y))/(max(y)- min(y))

plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
axis(side=2, labels=F)
rm(k)

fit.gnls<-gnls(y~b0+ ifelse(x- k <=0, b1*x, b1*x + b2*(x-k)),  start=c(k=15, b0=0.122, b1=0, b2 =0.05)) 

sig = summary(fit.gnls)$sigma # standard deviation of the residuals (called standard error in R)
lines(x, predict(fit.gnls))
lines(x, predict(fit.gnls)+sig, lty=2)
lines(x, predict(fit.gnls)-sig, lty=2) 

text(2, 1.45, "D", cex=1.2 )
text(3, 1.45, "Segmented linear relationship", pos=4, cex=0.8)
text(3, 1.25, "Constant variance", pos=4, cex=0.8)




# --------------------------------------------------------------
# Step increase in mean.
# --------------------------------------------------------------
x= seq(1,30, 1)
k= 15; d=15
y = ifelse(x-k<= 0, 10, 10 + d) + rnorm(length(x), 0, 5) 
y = (y-min(y))/(max(y)- min(y))

plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
axis(side=2, labels=F)
rm(k)
# rm(k) couldn't get model to converge on here. If I'm reading the paper correctly, need to give different values of k, calculate the likelihood and and then select the best model.
k1 = c(14, 15, 16)

logLikResult = rep(NA, 3)
for(i in 1:3){
  k = rep(k1[i], length(x))
  fit.gnls<-gnls(y~ ifelse(x- k <= 0, b0, b0 + b1),  start=c(b0=10,b1=15))
  logLikResult[i] = logLik(fit.gnls)
}
# k = 15 is best logLik score

k = rep(15, length(x))
fit.gnls<-gnls(y~ ifelse(x- k <= 0, b0, b0 + b1),  start=c(b0=10,b1=15))
sig = summary(fit.gnls)$sigma # standard deviation of the residuals (called standard error in R)
lines(x, predict(fit.gnls))
lines(x, predict(fit.gnls)+sig, lty=2)
lines(x, predict(fit.gnls)-sig, lty=2) 

text(2, 1.45, "E", cex=1.2 )
text(3, 1.45, "Step change relationship", pos=4, cex=0.8)
text(3, 1.25, "Constant variance", pos=4, cex=0.8)



# --------------------------------------------------------------
# Sigmoidal model (smooth step)
# --------------------------------------------------------------

M = 15 # mid point
B= 0.5 # slope
K = 40 # max
A = 0 # min
x = seq(1, 30, 1)

y = A+((K-A)/(1+ exp(-1*B*(x-M)))) +rnorm(length(x), 0, 5)
y = (y-min(y))/(max(y)- min(y))
fit.gnls<-gnls(y~A+((K-A)/(1+ exp(-1*B*(x-M)))), start=c(A=0.1 , K=0.9 , B=2, M=15))

plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
axis(side=2, labels=F)
sig = summary(fit.gnls)$sigma # standard deviation of the residuals (called standard error in R)
lines(x, predict(fit.gnls))
lines(x, predict(fit.gnls)+sig, lty=2)
lines(x, predict(fit.gnls)-sig, lty=2) 

text(2, 1.45, "E", cex=1.2 )
text(3, 1.45, "Smooth step change (sigmoidal)", pos=4, cex=0.8)
text(3, 1.25, "Constant variance", pos=4, cex=0.8)







### The challenges 

# --------------------------------------------------------------
# Linearly inc mean with linearly increasing variance.
# --------------------------------------------------------------
# 
# x= seq(1,30, 1)
# y = 10 + 2*x+ (x/10)*rnorm(length(x), 0, 5) 
# y = (y-min(y))/(max(y)- min(y))
# 
# plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
# axis(side=2, labels=F)
# fit.gnls<-gnls(y~b0+b1*(x),start=c(b0=0,b1=0) , weights =varFixed(~x))
# sig = attr(residuals(fit.gnls), "std")
# lines(x, predict(fit.gnls))
# lines(x, predict(fit.gnls)+sig, lty=2)
# lines(x, predict(fit.gnls)-sig, lty=2) 
# 
# text(2, 1.45, "C", cex=1.2 )
# text(3, 1.45, "Linear relationship", pos=4, cex=0.8)
# text(3, 1.25, "Increasing std.error", pos=4, cex=0.8)
# 






# --------------------------------------------------------------
# Step increase in mean and  variance
# --------------------------------------------------------------
# 
# x= seq(1,30, 1)
# k= 15; d=15
# y = ifelse(x-k<= 0, 10+ rnorm(length(x), 0, 5) , 10 + d + rnorm(length(x), 0, 15)) 
# y = (y-min(y))/(max(y)- min(y))
# 
# plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
# axis(side=2, labels=F)
# rm(k)
# 
# k1 = c(14, 15, 16)
# 
# logLikResult = rep(NA, 3)
# for(i in 1:3){
# k = rep(k1[i], length(x))
# fit.gnls<-gnls(y~ ifelse(x- k <= 0, b0, b0 + b1),  start=c(b0=10,b1=15))
# logLikResult[i] = logLik(fit.gnls)
# }
# # k = 15 is best logLik score
# 
# k = rep(15, length(x))
# 
# # split data into two sections
# varAB = rep(NA, length(x))
# varAB[which(x <= k[1])] = 1
# varAB[which(x > k[1])] = 2
# 
# fit.gnls<-gnls(y~ ifelse(x- k <= 0, b0, b0 + b1),  start=c(b0=10,b1=15), weights= varIdent(form=~1| varAB))
# sig = attr(residuals(fit.gnls), "std")
# lines(x, predict(fit.gnls))
# lines(x, predict(fit.gnls)+sig, lty=2)
# lines(x, predict(fit.gnls)-sig, lty=2) 
# 
# text(2, 1.45, "I", cex=1.2 )
# text(3, 1.45, "Step change relationship", pos=4, cex=0.8)
# text(3, 1.25, "Step change std.error", pos=4, cex=0.8)


#dev.off()
