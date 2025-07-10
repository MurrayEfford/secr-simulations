## Setup code for all secr-simulations

# Each set of simulations is in a folder secr-simulations/xxx
# where xxx is an alphanumeric code ('ARR', 'CLO' etc.) and the markdown file
# in that folder is 'secr-simulations-xxx.rmd'.

# This code is sourced by a startup chunk in each R Markdown file 

## To use
## -- install packages secrdesign and RColorBrewer from CRAN
## -- other packages are required for particular secr-simulations-xxx.rmd
## --     LDF  : mvtnorm
## --     OU   : mvtnorm
## --     RSF  : mvtnorm, RandomFields
## --     SARE : mvtnorm, raster, RandomFields
## -- edit the call to setNumThreads as required
## -- in the R Markdown (.rmd) file:
## --     edit nrepl
## --     set switches:
## --         runsim  if TRUE runs all simulation chunks
## --         figures if TRUE run figure chunks, assuming previous or current runsim
## --         cache   if TRUE cache selected chunks (mostly simulations)
## -- knit R Markdown file

library(secrdesign)
nc <- setNumThreads(min(parallel::detectCores(), 40))
options(digits = 5)      # for more readable output
options(width = 100)

# GENERAL

metadata <- function(runsim, extra = FALSE) {
    # runsim TRUE : novel simulations so we compute and store the simulation conditions
    # runsim FALSE: merely retrieve and include the previously stored conditions
    if (runsim) {
        # expect that internal code in calling function has set nrepl and possibly nrepl2
        nr1 <- if (exists("nrepl"))  paste0('nrepl      ', nrepl, '\n') else ''
        nr2 <- if (exists("nrepl2")) paste0('nrepl2     ', nrepl2, '\n') else ''
        md <- paste0( R.version.string, '\n',
                      'Platform   ', sessionInfo()$platform, '\n',
                      'Running    ', sessionInfo()$running, '\n',
                      'secr       ', packageVersion('secr'), ' ',  packageDate('secr'), '\n',
                      'secrdesign ', packageVersion('secrdesign'), ' ',  packageDate('secr'), '\n',
                      'Threads    ', setNumThreads(), '\n',
                      nr1, nr2,
                      'Run        ', date(), '\n')
        # optional supplements 2025-03-28
        if (extra) {
            md0 <- readRDS('metadata.RDS')
            md <- paste0(md0, 'Supplements\n', md)
        }
        saveRDS(md, file = 'metadata.RDS')
    }
    else {
        md <- readRDS('metadata.RDS')
    }
    cat(md)
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

# GRAPHICS
# ========

grdevice <- 'ragg_png'  # added 2025-07-10; used on NeSI

# predefine some colours
yob5 <- RColorBrewer::brewer.pal(n = 5, name = "YlOrBr")
blu5 <- RColorBrewer::brewer.pal(n = 5, name = "Blues")

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