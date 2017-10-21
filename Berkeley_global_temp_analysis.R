
## ======================================================================
## Global temperature series analysis
##
## "The Berkeley Earth project [Rohde et al., 2013] uses isotropic
## geostatistical tools (kriging) to estimate the global mean
## temperature. The land data have been collected from 14 databases
## and almost 45,000 stations. One of the main diﬀerences between the
## Berkeley Earth approach and most other global approaches is that
## the former group does not attempt to “homogenize” stations
## [Trewin, 2010]."
##
## Peter Craigmile and Peter Guttorp
## "Modeling and assessing climatic trends"
##
## Contact: pfc@stat.osu.edu
##
## Data source: A previous version of Land + Ocean (1850 – Recent),
## annual summary, downloaded in 2015 from
##   http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_summary.txt
## ======================================================================


## Load some useful R functions
source("functions.R")


## Read in the data
berkeley.temp <- read.table("Berkeley_Land_and_Ocean_summary.txt")

## Extract the necessary variables
the.mean <- berkeley.temp[,2] 
years    <- berkeley.temp[,1]
the.sd   <- berkeley.temp[,3]

the.vars <- the.sd^2


## Calculate the Bonferroni multiplier
Bonf <- abs(qnorm(0.025/length(the.mean)))




## Produce Figure 2

par(mfrow=c(1,2))

par(fig=c(0, 0.6, 0, 1))

anomaly.plot(ylim=c(-1.2, 0.8))

time.series.and.bands(years, the.mean, the.sd, multiplier=Bonf)

mtext("(a)", side=3, cex=0.8, line=0)

lines(years, the.mean)

par(fig=c(0.6, 1, 0, 1), new=TRUE)

plot(years, the.sd,
     xlab="Year",
     ylab="Estimated standard error (°C)", type="l")

mtext("(b)", side=3, cex=0.8, line=0)




## First load the nlme R library
library(nlme)

## Fit the OLS model, and summarize it
par(mfrow=c(1,1))

linear.ols <- gls(the.mean ~ years, method="ML")
summ("OLS", linear.ols)





## Produce Figure 5

par(mfrow=c(1,2))

par(fig=c(0, 0.6, 0, 1))

anomaly.plot(ylim=c(-1.2, 0.8))

time.series.and.bands(years, the.mean, the.sd, multiplier=Bonf)

lines(years, fitted(linear.ols))

mtext("(a)", side=3, cex=0.8, line=0)

lines(years, the.mean)

par(fig=c(0.6, 1, 0, 1), new=TRUE)

plot(years, resid(linear.ols),
     xlab="Year",
     ylab="OLS Residual (°C)", type="l")

abline(h=0, lty=2)

mtext("(b)", side=3, cex=0.8, line=0)



## Using orthogonal polynomials up to degree 2, define the design matrix
X <- cbind(1, poly(years, 2))
x1 <- X[,2]
x2 <- X[,3]



## Fit the weighted ARMA(4,1) model with the quadratic trend
weighted.ARMA41.quadratic <- gls(the.mean ~ x1 + x2, method="ML",
                                 correlation = corARMA(p=4, q=1),
                                 weights=varFixed(~the.vars))

summary(weighted.ARMA41.quadratic)




## install.packages("bentcableAR")

require(bentcableAR)

ts <- 1:length(the.mean)

## Fit the broken stick model

berk.bsar <- bentcable.ar(the.mean, p=4, stick=TRUE,
                          init.cable=c(0.4115256558,-0.0004087094, 0.0087471315,62))

berk.bsar.fit <- sapply(ts, function (t)
                        fullcable.t(t, berk.bsar$cable$est[1],
                                    berk.bsar$cable$est[2],
                                    berk.bsar$cable$est[3],
                                    berk.bsar$cable$est[4], 0))



## Fit the bent cable model

berk.bcar <- bentcable.ar(the.mean, p=4,
                          init.cable=c(-1.2677, 0.0004864, 0.0135112, 99, 66.4))

berk.bcar.fit <- sapply(ts, function (t)
    fullcable.t(t, berk.bcar$cable$est[1],
                berk.bcar$cable$est[2],
                berk.bcar$cable$est[3],
                berk.bcar$cable$est[4],
                berk.bcar$cable$est[5]))




## Produce Figure 6


## The estimated trend
fit <- predict(weighted.ARMA41.quadratic)

## The covariance matrix for the time series errors
V <- vcov(weighted.ARMA41.quadratic)

## Calculate the estimated standard error
se.fit <- sqrt(diag( X %*% V %*% t(X) ))

## and the Scheff\'e mutliplier
Scheffe <- sqrt(2* qf(0.95, 2, length(the.mean)-2))



par(mfrow=c(1,1))

anomaly.plot(ylim=c(-1.2, 0.8))

time.series.and.bands(years, fitted(weighted.ARMA41.quadratic), se.fit, Scheffe)

lines(years, the.mean)

lines(years, fitted(weighted.ARMA41.quadratic))

lines(years, berk.bsar.fit, lty=2)

lines(years, berk.bcar.fit, lty=3)

legend(1970, -0.45, c("Quadratic", "Broken Stick", "Bent Cable"),
       lty=c(1,2,3), cex=0.85, bty="n")






library(splines)

## Fit the weighted ARMA(2,1) model with the df=8 spline
weighted.ARMA21.df8.spline <- gls(the.mean ~ bs(years, 8), method="ML",
                                 correlation = corARMA(p=2, q=1),
                                 weights=varFixed(~the.vars))

summary(weighted.ARMA21.df8.spline)


## Fit the weighted ARMA(4,1) model with the df=5 spline
weighted.ARMA41.df5.spline <- gls(the.mean ~ bs(years, 5), method="ML",
                                 correlation = corARMA(p=4, q=1),
                                 weights=varFixed(~the.vars))

summary(weighted.ARMA41.df5.spline)



pdf(file="../trend_chapter_r1/figures/fig6_nonlinear_trends.pdf", width=4.6, height=2)
par(mfrow=c(1,1), cex=0.7, mar=c(2.7,2.7,1.3,0.2), mgp=c(1.7,0.5,0), bty="L")

anomaly.plot(ylim=c(-1.2, 0.8))

time.series.and.bands(years, fitted(weighted.ARMA41.quadratic), se.fit, Scheffe)

lines(years, the.mean)

lines(years, fitted(weighted.ARMA41.quadratic), lwd=1.5)

lines(years, fitted(weighted.ARMA41.df5.spline), lty=2, lwd=1.5)

lines(years, fitted(weighted.ARMA21.df8.spline), lty=5, lwd=1.5)


legend(1950, -0.45, c("Quadratic + ARMA(4,1)",
                      "Spline, df=5 + ARMA(4,1)",
                      "Spline, df=8 + ARMA(2,1)"),
       lty=c(1,2,5), cex=0.85, bty="n", lwd=1.5)

dev.off()
