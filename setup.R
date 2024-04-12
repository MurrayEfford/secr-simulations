## Setup code for all secr-simulations
## 2024-03-29

## To use
## -- install package secrdesign from CRAN
## -- edit the call to setNumThreads as required

library(secrdesign)
nc <- setNumThreads(18)
options(digits = 5)      # for more readable output
options(width=100)

# add points to plot
addRB <- function (x, est, mean = 'RB', se = 'seRB', xoffset = 0, star = 100, ...) {
    x <- x+xoffset
    OK <- est[,mean] <= star
    OK[is.na(OK)] <- FALSE
    points(x[!OK], rep(star,sum(!OK)), pch=8)
    mn <-  est[OK,mean]
    se <-  est[OK,se]
    segments(x[OK],mn-2*se, x[OK], mn+2*se)
    points(x[OK], mn, ...)
}

# add shaded strip for +/- RB%
shade <- function(RB = 0.1) {
    x <- par()$usr[1:2]
    if (par()$xlog) x <- 10^x
    polygon(x = c(x,rev(x)), y= RB * c(-1,-1,1,1), col= grey(0.92), border = NA)
    box()
}

# recent secrdesign automatically saves summary capture statistics when fit = TRUE
# this function retrieves a summary for each scenario
getcounts <- function (sims, fn = mean) {
    test <- attr(sims$output[[1]][[1]], 'counts')
    if (is.null(test))
        stop ("counts were not saved with these simulations - only available for secrdesign>=2.9.1")
    summarycounts <- function(x) {
        countmatrix <- sapply(x, attr, 'counts')
        if (is.matrix(countmatrix))
            apply(countmatrix,1,fn)
        else
            rep(NA,length(test))
    }
    t(sapply(sims$output, summarycounts))
}

countlegend <- function (rows = -5) {
    df <- data.frame(Variable = c('n','r','nmov','dpa','rse','rpsv'), 
                     Definition = c(
                         'Number of individuals', 
                         'Number of recaptures', 
                         'Number of movements', 
                         'Detectors per individual',
                         'Approximate RSE',
                         'Approximate sigma-hat'))
    print(df[rows,], row.names = FALSE, right = FALSE)
}
