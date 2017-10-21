
## ======================================================================
## Useful R functions for:
##
## Peter Craigmile and Peter Guttorp
## "Modeling and assessing climatic trends"
##
## Contact: pfc@stat.osu.edu
## =======================================================================

time.series.and.bands <- function (x, y, the.sd,
                                   multiplier=1.96, bands=TRUE) {

    ## Produces a line plot of y versus x along with confidence bands or lines.

    pm <- multiplier*the.sd

    if (bands) {
        
        polygon(c(x, rev(x)), c(y-pm, rev(y+pm)), col="gray70", border=NA)
    } else {
        
        segments(x, y-pm, x, y+pm, col="gray70", lwd=2) 
        
    }

    lines(x, the.mean, lwd=1)
}



anomaly.plot <- function (ylim=c(-1.2, 1.2)) {
    
    plot(years, the.mean, 
         xlab="Year",
         ylab="Temperature anomaly (°C)", type="n",
         ylim=ylim, yaxt="n")
    
    axis(side=2, at=seq(-1, 1, 0.5))
}


AIC.to.AICC <- function (aic, n, npars) {
  ## ======================================================================
  ## Purpose : Transforms the AIC value to the AICC value.
  ## ======================================================================
  
  aic - 2 * npars * ( 1 - n/(n-1-npars))
}




summ <- function (the.name, model) {
    ## ======================================================================
    ## Summarize the 'model' giving it the name 'the.name'
    ## Produces a table in LaTeX format; i.e., with "&" between each column.
    ## ======================================================================
        
    res <- resid(model, type="n")

    t.table <- (summary(model))$tTable

    beta1.hat    <- t.table[2,1]*100
    se.beta1     <- t.table[2,2]*100
    beta1.Pvalue <- t.table[2,4]

    if (is.character(names(model$modelStruct))) {

        npars <- 0
        
    } else {

        npars <- length(AR1$modelStruct)
    }

    the.AIC      <- AIC(model)

    #AIC.to.AICC(AIC(model), model$dims$N,
    #                            model$dims$p+1+npars)
    
    Lj.Pvalue    <- Box.test(res, 10, type="Lj")$p.value

    acf(res, main="")

    cat(the.name)
    cat(" & ")
    cat(sprintf("%1.2f", beta1.hat))
    cat(" & ")
    cat(sprintf("%1.2f", se.beta1))
    cat(" & ")
    val <- sprintf("%1.4f", beta1.Pvalue)
    if (val=="0.0000") {
        cat("$<$0.0001")
    } else {
        cat(sprintf("%1.4f", beta1.Pvalue))
    }
    ##    cat(signif(beta1.Pvalue, digits=2))
    cat(" & ")
    cat(sprintf("%1.1f", the.AIC))
    cat(" & ")
    cat(sprintf("%1.3f", Lj.Pvalue))
    cat(" \\\\")
    cat("\n")
}

