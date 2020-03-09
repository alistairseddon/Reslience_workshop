####
library("nlme")
# One of the classic pieces of evidence used to assess regime shifts in palaeoecological temporal series are abrupt changes in the response variable (see Scheffer and Carpenter 2003)

# There a number of ways to do this:

# Classical methods (e.g. constrained cluster analysis) and temporal series type approaches.
# Here we focus on some commonly used approaches assessing temporal series
# - We can look at Sequential T-Test analysis (STARS, Rodionov et al. 2003,2006)
# - Break-point analysis using non-linear regression (Carstensen and Weydmanm 2011, used in Seddon et al. 2014)
# - Generalised additive modelling (Simpson 2019)

## Simulated data to test and explore the approaches

# Simulate some data with a clear shift in mean
t <- 1:100
A <- rnorm(length(t))
A[51:100] <- A[51:100] + 1
A[1:50] <- A[1:50] - 1
simulatedData <- tibble (t = t, A = A)

regime1_mean <- round(mean(A[1:50]), 2)
regime2_mean <- round(mean(A[51:100]), 2)

fit.gnls<-nlme::gnls(A~ ifelse(t- 50 <= 0, b0, b0 + b1),  start=c(b0=-2,b1=2))

# Check this value against the regime1_mean/ regime2
fit.gnls
predicted.values <- predict(fit.gnls)

# Here we know there is a clear shift in the mean value
simulatedData %>% 
  ggplot(aes(x = t, y = A)) +
  geom_point() +
  geom_line(y = predicted.values, color = "red") +
  annotate(geom = "text", x=25, y=0, label=paste0("R1_mean = ", regime1_mean),  color = "blue") +
  annotate(geom = "text", x=75, y=-0, label=paste0("R2_mean = ", regime2_mean), color = "blue")
  
# So what are some ways to test this?

# STARS was an algorithm proposed by Rodionov (2003)

# Take figure from Rodionov paper

PDO <- read_csv("data/PDO.csv")
PDO_new <- t(sapply(PDO$Date, function(x) substring(x, first=c(1,5), last=c(4,6))))
colnames(PDO_new) <- c("Year", "month")
PDO <-  cbind(PDO_new, PDO[,"Value"])
write.csv(x = PDO, file = "data/PDO_NOAA.csv", row.names = FALSE  )

# PDO <- read_csv("data/PDO_NOAA.csv") %>% as.tibble() %>% filter(month == "01") %>% filter(Year >= 1900)

PDO <- read.delim("data/PDO_other.txt", sep = " ") %>% 
  as.tibble() %>% 
  filter(YEAR <= 2003)
View(PDO)

PDO_plot <- PDO %>% 
  select(YEAR, JAN) %>% 
  filter(YEAR <= 1927) %>% 
  ggplot(aes(x = YEAR, y = JAN)) +
  geom_bar(stat = "identity")

# take a given window length : l, and a given p value (e.g. 0.05)
#  1. Set the cut off length l of the regimes for variable, A  (more information on L given below)

l <- 10
A <- PDO$JAN
#  2. find the difference (diff) between mean values of two subsequnt regimes that would be statistically signficant according to students' t-test
# diff = t * sqrt(2 sigma<l>^2)/ l
# t= value of t distribution with 2l  2 degrees of freedom, assuming equal variance for both regimes and same for the average variance of the running l-year intervals in the time series of variable X

# moving window to estimate running variance of January PDO
var.l <- rep(NA, (length(A)-l+1)) # Set an empty results vector filled with NA
for(i in 1: length(var.l)) { # initialise the loop
  A_subset <- A[i : (i + l - 1)] # subset the window that you want
  var.l[i] <- var(A_subset) # get the variance
}
runningVariance <- mean(var.l) # This is a slightly different number to that reported in Rodionov?
runningVariance
runningVariance <-  0.76 # Set it to the value he reported in paper
t <- qt(p = 0.025, df = 2*l-2) *-1 # Two tailed
# diff = t* sqrt(2 * runningVariance / l)
diff = t* sqrt(2 * runningVariance / l)

# 3. Calculate the mean (A) r1 of the initial values of A. Levels that should be reached in subsequent l years to qualify for shift to regime R2 mean(A)r2 = mean(A)r1± diff
# Mean value for first 10 year window
r1.mean <- mean(A[1:l])
r1.mean

# Mean value for regime R2, should be greater than r1.mean + diff or less than this
r2.mean.level_upper <- r1.mean + diff # should be greater than this value (1.42709)
r2.mean.level_lower <- r1.mean - diff # or less than this value (-0.2110902)
r2.mean.level_upper
r2.mean.level_lower

# 4. For each value of starting with year i = l +1, see if it is greater than or less than the critical value
# E.g in 1910 the value of the Jan PDO was -2.33
PDO %>% filter(YEAR == 1910)
# Value was -0.025

# NB, values are the same signs so Rodionov multplies by *-1 to get the final RSI here. Only do this if going from negative to positive
RSI_1910_1912 <- ((-0.25 - r2.mean.level_lower) / sqrt(runningVariance)/ l )*-1
# Can turn this into a function
RSI.calc.lower <- function(value, threshold){
  # value = next value to test in window for regime shift
    RSI <- ((value - threshold) / sqrt(runningVariance)/ l )*-1
 return(RSI)
}

RSI.calc(-0.25)

PDO %>% filter(YEAR == 1911) # -1.11
PDO %>% filter(YEAR == 1912) # -1.72

RSI_1910_1912 <-  RSI.calc.lower(-0.25) + RSI.calc.lower(-1.11) +  RSI.calc.lower(-1.72) 
RSI_1910_1912 

# Or you could apply it to the whole window
i <- 11
PDO_window <- PDO[i: (i+ l-1), c("YEAR","JAN")]

PDO_window <- PDO_window %>% 
  mutate(RSI = map2(.x = JAN, .y = r2.mean.level_lower, .f = RSI.calc.lower)) %>% 
  unnest(cols = c(RSI)) %>% 
  mutate(cumsum=cumsum(RSI))
  
PDO_window %>% 
  ggplot(aes(x = YEAR, y = cumsum)) +
  geom_bar(stat = "identity")

# Find the whole RSI value. See it is positive so we confirm the regime shift at this starting value
sum(PDO_window$RSI)


# Search for a new regime shift, to (R3) starts at 1911 (i = 12), using the mean value of those first 10 samples in the regime R2 as the mean value
r2.mean <- mean(PDO_window$JAN)
r3.mean.level_upper <- r2.mean + diff
r3.mean.level_lower <- r2.mean - diff

i <- 12
newValueToTest <- PDO[i, "JAN"]
newValueToTest > r3.mean.level_upper
newValueToTest < r3.mean.level_lower

# Therefore RSI = 0, we move on

i <- 13
newValueToTest <- PDO[i, "JAN"]
newValueToTest > r3.mean.level_upper
newValueToTest < r3.mean.level_lower

PDO_window <- PDO[i: (i+ l-1), c("YEAR","JAN")]

PDO_window <- PDO_window %>% 
  mutate(RSI = map2(.x = JAN, .y = r3.mean.level_lower, .f = RSI.calc.lower)) %>% 
  unnest(cols = c(RSI)) %>% 
  mutate(cumsum=cumsum(RSI))

PDO_window %>% 
  ggplot(aes(x = YEAR, y = cumsum)) +
  geom_bar(stat = "identity")
# RSI becomes negative after one year so we move on.

# According to the paper we don't see another regime shift until year 1922 when the new regime is identified

# If not true, move onto next value in the dataset and retest

# If true, each new value of xi, whee i > j, is used to confirm or reject the null hyphtesis of a regime shift at year j.
# In this case you calculate the regime shift index (RSI), which is the cumulative sum of the normalised anomalies:

# Do this for every window- if at any time this turns negative, proceed to 6. If positive, move to 7

# Step 6: Negative RSI value means test for a regime shift at year j failed. Assign RSi to 0 and start again. Recaluate xR1 to inlcude the value of xj and keep testing the values of xi starting with i = j+1 

# Step 7. Positive value of rSI means that the regime shift at year j is significant at the probability level p. Calculate the actual mean value for R2. In the first instance it is just one value, but this conitues to find the next regime.
# Calculations continue i a loop until all data for X are processed.

# Thankfully we already have some code that can be used that does all of this.

source("STARS_update.R")
# Run stars
PDO.RS <- stars(y= PDO$JAN, time  =  PDO$YEAR, L=10, p=0.1, h =1,  AR1red="none", prewhitening = FALSE)
# Potentially a bug in the coding which means 1 tailed not two tailed test is used? Need to check this

# Plot result
plot.stars(PDO.RS) # NB, gives slightly different results to the Rodionov algoirthm, I think because of the slightly different weighting functions.

# Can inspect the object produced by the function here
PDO.RS

## Highly dependent on window length and p value
# If window length increases, the number of degrees of freedom increases and so the statistically significant difference needed between two regimes becomes smaller
# Translates as higher values of RSI for same  value. 
# Window size influences detection capabilty of different regimes of different magnitudes


# Gaussian, white noise process. Does not take autocorrelation into account.


### Autocorrelation section here

# Rodionov builds in some correlation functions here to do this

PDO.RS_ar <- stars(y= PDO$JAN, time  =  PDO$YEAR, L=10, p=0.1, h =1,  AR1red="IP4", prewhitening = TRUE)
plot.stars(PDO.RS_ar) # Doesn't work needs fixing



#### Problems for palaeoecological data?
#- needs evenly spaced time series
#- very dependent on moving window size
#- What about alternatives?

### Generalised non-linear regression is a useful tool for analysing break-points in time series

# Fitting a model by ordinary least squares
# aims to miminise the difference between predicted value and the residuals.
# Makes certain assumptions:
# 1. The relationship between the response and the predictors is ~linear.
# 2. The residuals have a mean of zero.
# 3. The residuals have constant variance (not heteroscedastic).
# 4. The residuals are independent (uncorrelated).
# 5. The residuals are normally distributed.

# Example

b1 <- 2
t <- 1:20
y <- b1* t + rnorm(length(t), 0, 5)

ols <- lm(y ~ t)
pred_values <- predict(ols)

plot(t, y, ylim = c(0,50))
abline(ols, col = "blue")
ols.data <- cbind(t, y, pred_values)
for(i in 1: nrow(ols.data))  lines(c( ols.data[i,"t"], ols.data[i,"t"]), c(ols.data[i,"y"], ols.data[i,"pred_values"]), col= "red")
residuals <- y - predict(ols)
hist(residuals)

coef(ols)
plot(ols)
vcov(ols)

## But, with time series data, as we have seen, residuals may not be autocorrelated. They may also not have constant variance. What other approaches are there?

# One possible solution; Generalised (weighted) regression GLS (generalised least squares)

# y.i = b0 + b1*x.i + ε.i

# Σ ε.i ^2 = Σ { y.i - (b0 + b1*x1)} ^2

# y - response
# x = predictor

# W weight matrix

weight.matrix <- matrix(0, ncol = length(t), nrow = length(t))
diag(weight.matrix) <- 1

# We could downweight an observation because of sample error 

# e.g. downwieght 5th observation

weight.matrix[5,5] <- 0.5
weight.matrix
###


# Why is this useful, because residuals in a dataset might not necessarily have uniform variance

# Example. Now imagine we have the same model as before, but the size of the error increases with the size of t
b1 <- 2
t <- 1:20
y <- b1* t + rnorm(length(t), 0, 5)*t/10

ols <- lm(y ~ t)
pred_values <- predict(ols)

plot(t, y, ylim = c(0,50))
abline(ols, col = "blue")
ols.data <- cbind(t, y, pred_values)
for(i in 1: nrow(ols.data))  lines(c( ols.data[i,"t"], ols.data[i,"t"]), c(ols.data[i,"y"], ols.data[i,"pred_values"]), col= "red")
residuals <- y - predict(ols)

plot(t, residuals)
abline(h = 0, lty = 2)

# So our ols fit is biased
summary(ols)

# We can use gls to weight the residuals in the model
#kiyi = kiβ0 + kiβ1xi + kiεi
# W = weight matrix
# k = sqrt(K) #is it always sqrt(W)??
# W = kk

# So every parameter value is multiplied by the apppropriate weight. We can therefore adjust the influence of the point and also adjust the szie of the residuals

# Can be estimated using gls funciton in R

gls_fit <- gls(y ~t )
gls_fit_weights <- gls(y ~t, weights =varFixed(~t) )

summary(gls_fit)
summary(gls_fit_weights)

# Can check the weights here
modelWeights <- attr(gls_fit_weights$modelStruct$varStruct, "weights")

# Can get the weighted residuals 
w_resid <- resid(gls_fit_weights, "normalized") # gls_fit_weights$residuals*(modelWeights)
un_weighted_resid <- resid(gls_fit_weights) 

par(mfrow = c(1,1))
plot(t, un_weighted_resid)
abline(h= 0, lty = 2)
points(t,  w_resid, col = "blue" )  # Looks much better

# We can use AIC to check whether this is in fact a better model
AIC(gls_fit, gls_fit_weights ) # The model with the weighted residuals is better due to lower AIC

# Why are there the same numbers of degrees of freedom here?

# There are many different variance structures to explore
gls_fit_weights_2 <- gls(y ~t, weights =varExp(form = ~t) )
AIC(gls_fit, gls_fit_weights, gls_fit_weights_2 ) 

w_resid_2 <- resid(gls_fit_weights_2, "normalized")
par(mfrow = c(1,1))
plot(t, un_weighted_resid)
abline(h= 0, lty = 2)
points(t,  w_resid, col = "blue" )  # Looks much better
points(t,  w_resid_2, col = "red" )  # More complicated model, don't get much difference in error structure


### Also useful for, e.g. populations with different variances
# e.g. could be useful for, e.g., test for changing variance over time in a temporal series

r1 <- rnorm(10, 0, 0.5)
r2 <- rnorm(10, 0, 3)
y <- c(r1, r2)
regime <- as.factor(c(rep("r1", 10), rep("r2", 10)))
plot(t, y,col = regime)
 
gls_unweighted <- gls(y ~ t)
gls_weighted <- gls(y ~ t, weights = varIdent(form = ~+1|regime))

AIC(gls_unweighted, gls_weighted)
summary(gls_unweighted)

#### Temporal autocorrelation.

# gls can also be used when modelling temporal autocorrelation

K = 100  # Carrying capacity, range 50-100
r = 0.2 # Growth rate, range 0-1

# Consumer parameters
g = 4 # maximum consumption rate
H = 5 # Sat. rate determines centre point of sigmoidal curve
p = 5  # power term, gives sigmoidal relationship and determines steepness
#c =  1   # number of consumer


  t <- 1: 180
  A <- rep(NA, length(t))
 # K <- 10 # Carrying capacity
  # cmax <- 3.9 # Max number of consumer (1-3.9)
  #crange <- 1 #range of c(0-2)
  cmax <- c( seq(0.1, 2, 0.02), rep(2, 79))
    
  
  starting.A <- 90
  
  A[1] <- starting.A
  dA.dt <- rep(NA, length(t)) # Rate of population growth for timestep 1
  
  # for(i in 1: (length(t)-1)) {
  #   dA.dt[i] <-   A[i] * (1 -(A[i]/K) )  - runif(min = 0.3, max =0.8, n =1)*cmax[i]*(A[i]^2 / (A[i]^2 + 1))
  #   A[i+1] <- A[i]+ dA.dt[i] + rnorm(1, 0, 0.8)
  # }
  # 
  # crash.model <- tibble(t = t, A= A, At.1 = c(A[-1], NA))[1:180,]
  # 
  for(i in 1: (length(t)-1)) {
    dA.dt[i] <-  r* A[i] * (1 -(A[i]/K) )  -  runif(min = 0.1, max =0.9, n =1)*(cmax[i]*g *A[i]^p / (A[i]^p + H^p))
    A[i+1] <- A[i]+ dA.dt[i]
  }  
  crash.model2 <- tibble(t = t, A= A, At.1 = c(A[-1], NA))# [45:180,]
  
 plot(crash.model2$t, crash.model2$A)

save(crash.model2, file ="crash.model2.RData")
load("crash.model2.RData")

t <- 1:132
y <- crash.model2$A[1:132]
regime <- as.factor(c(rep("r1", 33), rep("r2", 99)))
plot(t, y, col = regime)

# test if there is a regime shift at time t = 36

gls_null <- gls(  y ~ 1) # null model

gls_mod_1<-gls(y ~ regime ) # model between the different components
summary(gls_mod_1)
# yes, thinks there is a significant difference in means between these two time points

# but check the residuals
plot(t, resid(gls_mod_1)) # definitely looks like there is some autocorrelation in the structure of the model

# Use the acf function to test if there is autocorrelation in the residuals of the model. Yes there is. So our model might be overfit
# what's the difference between partial on full lag one autocorrelation
acf(resid(gls_mod_1)) # Indeed there is

# OK so now we can fit a model which takes this into account.

# Can use the nlme corAR1 function to respecify the correlation matrix which takes these autocorrelated errors into account
gls_mod_2<-gls(y ~ regime,  correlation = corAR1(form = ~ t))
summary(gls_mod_2)

gls_mod_2_resid <-  resid(gls_mod_2, type ="normalized")
acf(gls_mod_2_resid) # Much better

# Ok but is there actually a hint of increase in variance as well?
plot(t, gls_mod_2_resid) 

gls_mod_3<-gls(y ~ regime,  correlation = corAR1(form = ~ t), weights = varIdent(form = ~+1|regime))
summary(gls_mod_3)
ACF(gls_mod_3, resType = "normalized", form = ~ t)

gls_mod_4<-gls(y ~ regime,  correlation = corAR1(form = ~ t), weights = varFixed(~t))
summary(gls_mod_4)


gls_mod_5<-gls(y ~1 ,  correlation = corAR1(form = ~ t))
summary(gls_mod_5)

AIC( gls_mod_1, gls_mod_2, gls_mod_3, gls_mod_4, gls_mod_5) # Check error message

gls_mod_3_resid <-  resid(gls_mod_3, type ="normalized")
plot(t, gls_mod_3_resid) 

str(summary(gls_mod_3)$ modelStruct $ corStruct)

# What is it actually doing?

# The correlation matrix of a model which assumes indepedence between datapoints
cs0 <- corAR1(0)
corMatrix(cs0, y)

# The correlation matrix of a model with a high order of autocorrelation (phi = 0.98). The AR1 (first order auto-regressive structure)
# defines a correlation structure in whihc the degree of correlation between the two observations/ residuals is proportional to 
# the relative amount of elapsed time. Specifically, the coorelation is defined as phi^|t-s| wheere |t-s| is the absolute diference
# between the current time (t) and previous time (s). So correlation of lag one is Phi, lag 2 is phi^2 etc.
# This autocorrelation structure assumes that (i) the pattern of correlation decline depends 
# only on the lag and not when in the timeseries the lag occurs (i.e. stationarity)
# (ii) correlation decline follows a very specific exponential pattern

cs1 <- corAR1(0.98, form =1 ~ t)
corMatrix(cs1, y)

# Here, the value phi is the autocorrelation coefficient.

# https://www.flutterbys.com.au/stats/tut/tut8.3a.html
# So what gls is doing here is refitting the model with the different correlation structure in the model. 
# Using the correlation=corAR1() argument, we can indicate that the off-diagonals of the variance-covariance matrix
# should follow a specific pattern in which the degree of correlaton varies as an exponent equal to the lag. 
# The spacing of the lags is taken from the form= argument. Of course, conceptually this is much more straight forward. 
# If this is left blank, then it is assumed that the data are in chronological order and that each row represents a single unit increase in time (lag). 
# In this case, the year variable is just a vector containing integers from 1 to 50, and therefore not strictly required. 
# Nevertheless, it is good practice to always include it.


# Incidentally, what happens if we assessed our timeseries with STARS

crash.stars <- stars(y= y, time  =  t, L=20, p=0.05, h =1,  AR1red="none", prewhitening = FALSE)
crash.stars_ar <- stars(y= y, time  =  t, L=20, p=0.05, h =1,  AR1red="IP4", prewhitening = TRUE)
plot.stars(crash.stars) 



### One more thing. Non linear regression: (GNLS) gives you the opportunity to fit many different model forms

# The same principles as GLS except here we can specify different model functions
# Proposed in Carstensen and Weydmenn (2012) for modelling time series.
# Breakpoints, step changes in mean, step changes in variance, segemented changes in mean

# all can be tested.

# Approach is very similar, but use a slightly different model specification and algorthm in gnls

# Example # Linear model fit
# Fit the data

x <- seq(1,30, 1)
y <- 10 + (2*x) + rnorm(length(x), 0, 5) 
y <- (y-min(y))/(max(y)- min(y))

fit.gnls<-gnls(y~ b0 + b1*(x),start=c(b0=0, b1=0) ) # You have to write the equation yourself, and also provide starting values for the iteration to fit

# Make a plot
plot(x, y, pch=20, ylim=c(-0.5,1.5), main = "", yaxt="n", xaxt="n")
axis(side=2, labels=F)
sig <- summary(fit.gnls)$sigma # standard deviation of the residuals
lines(x, predict(fit.gnls))
lines(x, predict(fit.gnls)+sig*2, lty=2)
lines(x, predict(fit.gnls)-sig*2, lty=2) 

text(2, 1.45, "A", cex=1.2 )
text(3, 1.45, "Linear relationship", pos=4, cex=0.8)
text(3, 1.25, "Constant std.error", pos=4, cex=0.8)

# See the attached the R script. Multiple responses modelled to show different relationships
source("NLS_Script_ModelFitting.R")

# Can you modify the code to make a step change model with both a step change relationship and a step change in variance at the same time?
# Can you modify the code to show a linear model with an increase in variance?


### Now on to the palaeoecology. Exercise.

# Aims- to download data from Neotoma and fit a model to test for the temporal dynamics of pollen communities

# Step 1: Download XXX site from Neotoma and calculate an ordination on pollen data to describe the main trends in vegetation through time
# Step 2: Use gnls modelling to investigate the temporal dynamics across XX. Is there a trend in the vegetation composition?
# Step 3: What extra data processing techniques would be required to use the STARS algorithm?
# Step 4: How does this appraoch compare to other, classical modelling approaches?

# Distantia?















