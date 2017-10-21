
## ======================================================================
## US average temperature anomaly analysis in R
##
## "Shen et al. [2012] produced statistical estimates of US
## temperature anomalies from 1897--2008, using the US Historical
## Climatology Network data set version 2 [Menne et al., 2009],
## corrected for the fact that the time of day that measurements are
## made at can diﬀer by site."
##
## Peter Craigmile and Peter Guttorp
## "Modeling and assessing climatic trends"
##
## Contact: pfc@stat.osu.edu
##
## Data source:
##   Dr. Samuel Shen
##   Department of Mathematics and Statistics
##   San Diego State University 
## ======================================================================

## Load some useful R functions
source("functions.R")


## Read in the US average temperature anomalies
## First line is header
## No SD for 1895-6 so we don't use them
usyearly <- read.table("Fig6_ts_table1.txt", header=FALSE, skip=3)

## Create the necessary variables
years    <- usyearly[,1]
the.mean <- usyearly[,2]
the.sd   <- usyearly[,3]

## Calculate the variances
the.vars    <- the.sd^2



## Calculate the Bonferroni multiplier
Bonf <- abs(qnorm(0.025/length(the.mean)))



## Produce Figure 1

par(fig=c(0, 0.6, 0, 1))

anomaly.plot(ylim=c(-1.1, 1.1))

time.series.and.bands(years, the.mean, the.sd, multiplier=Bonf)

mtext("(a)", side=3, cex=0.75, line=0)

par(fig=c(0.6, 1, 0, 1), new=TRUE)

plot(years, the.sd,
     xlab="Year",
     ylab="Estimated standard error (°C)", type="l")

mtext("(b)", side=3, cex=0.75, line=0)



## Produce Table 1:

## First load the nlme R library
library(nlme)


## Fit the OLS model, and summarize it
## Format of summary is:
## Model name & Slope per century & Standard error & Slope P-Value & AIC & Ljung-Box P-value
linear.ols <- gls(the.mean ~ years, method="ML")
summ("OLS", linear.ols)


## WLS model
linear.wls <- gls(the.mean ~ years, method="ML",
                  weights=varFixed(~the.vars))
summ("WLS", linear.wls)


## AR(1) and weighted AR(1)
AR1 <- gls(the.mean ~ years, method="ML",
           correlation = corAR1())
summ("AR(1)", AR1)

weighted.AR1 <- gls(the.mean ~ years, method="ML",
                    correlation = corAR1(), weights=varFixed(~the.vars))
summ("Weighted AR(1)", weighted.AR1)


# AR(4) and weighted AR(4)
AR4 <- gls(the.mean ~ years, method="ML",
           correlation = corARMA(p=4, q=0))
summ("AR(4)", AR4)

weighted.AR4 <- gls(the.mean ~ years, method="ML",
                    correlation = corARMA(p=4, q=0), weights=varFixed(~the.vars))
summ("Weighted AR(4)", weighted.AR4)



## ARMA(3,1) and weighted ARMA(3,1)
ARMA31 <- gls(the.mean ~ years, method="ML",
              correlation = corARMA(p=3, q=1))
summ("ARMA(3,1)", ARMA31)

weighted.ARMA31 <- gls(the.mean ~ years, method="ML",
                       correlation = corARMA(p=3, q=1), weights=varFixed(~the.vars))
summ("Weighted ARMA(3,1)", weighted.ARMA31)





## Produce Figure 3

par(mfrow=c(1,2))

acf(residuals(linear.ols, type="n"), ylim=c(-0.2, 1), main="")
mtext("OLS residuals", side=3, line=0, cex=0.7)

acf(residuals(linear.wls, type="n"), ylim=c(-0.2, 1), main="")
mtext("WLS residuals", side=3, line=0, cex=0.7)





## Produce Figure 4

## The estimated trend
fit <- predict(weighted.ARMA31)

## The covariance matrix for the time series errors
V <- vcov(weighted.ARMA31)

## The design matrix
X <- cbind(1, years)

## Calculate the estimated standard error
se.fit <- sqrt(diag( X %*% V %*% t(X) ))

## and the Scheff\'e mutliplier
Scheffe <- sqrt(2 * qf(0.95, 2, length(the.mean)-2))


par(mfrow=c(1,1))

anomaly.plot(ylim=c(-1.1,1.1))

time.series.and.bands(years, fitted(weighted.ARMA31), se.fit, Scheffe)

lines(years, the.mean)

lines(years, fitted(linear.wls), lty=3, lwd=1)

lines(years, fitted(weighted.AR4), lty=2, lwd=1)

legend(1970, -0.45, c("WLS", "Weighted ARMA(3,1)", "Weighted AR(4)"),
       lty=c(3,1,2), cex=0.85, bty="n", lwd=1)
