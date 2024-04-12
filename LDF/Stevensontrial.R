## Derived from code of Ben Stevenson
## also generates trps, trps2, msk, msk2

convertCH <- function (capt, traps) {
    n <- nrow(capt)
    ch <- array(capt, dim=c(n,1,ncol(capt)), dimnames=list(1:n,NULL,NULL))
    class(ch) <- 'capthist'
    traps(ch) <- traps
    ch
}

sim.sscr.mge <- function(
        popn,
        traps, 
        cov.structure = "none", 
        det.pars = NULL, 
        cov.pars = NULL,
        output = c('pdot','umat','hazard','capt')
        ) {
    
    output <- match.arg(output)
    resp.pars <- 1
    detfn <- "hhn"
    re.multiplier <- "er"  ## model the hazard
    
    ## Sorting out identifiability.
    if (!(cov.structure %in% c("none", "OU")) & re.multiplier == "er" & detfn == "hhn"){
        if (is.null(cov.pars$mu.u)){
            cov.pars$mu.u <- 0
        } else {
            if (cov.pars$mu.u != 0){
                warning("The mu.u and lambda0 parameters are not identifiable for models with a hazard halfnormal detection function and an encounter-rate random effect multiplier. Setting nonzero mu.u is not recommended.")
            }
        }
        if (is.null(det.pars$lambda0)){
            stop("A value for lambda0 must be supplied.")
        }
    }
    
    ## Extracting number of traps.
    n.traps <- nrow(traps)
 
    ## number of individuals.
    n <- nrow(popn)
 
    ## Distances from activity centres to traps.
    ac.dists <- secr::edist(popn,traps)
   
    if (cov.structure == "OU"){
        if (output != 'capt') stop ("OU is only for capt simulation")
        ## Extracting movement parameters.
        tau <- cov.pars$tau
        sigma <- cov.pars$sigma
        n.steps <- cov.pars$n.steps   # effective "duration"
        epsilon <- cov.pars$epsilon   # proximity threshold for detection
        capt <- matrix(0, nrow = n, ncol = n.traps)
        for (i in 1:n){
            locs <- sim.ou(as.numeric(popn[i, ]), tau, sigma, n.steps)
            capt[i, ] <- count.dets(locs, traps, epsilon)
        }
        capt <- capt[apply(capt, 1, function(x) sum(x) > 0), ]
        convertCH (capt, traps)
    } else {
        
        ## Sorting out detection function.
        calc.detfn <- detfn.closure(detfn, det.pars)
        
        ## Distances between traps.
        trap.dists <- as.matrix(dist(traps))
        
        ## Calculating "baseline" encounter rates 
        base.er <- calc.detfn(ac.dists, er = TRUE)
        
        ## Constructing covariance matrix.
        if (cov.structure == "none"){
            mu.u <- 0
            cov <- matrix(0, nrow = n.traps, ncol = n.traps)
        } else if (cov.structure == "independent"){
            ## Extracting parameters.
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            ## Specifying covariance.
            cov <- matrix(0, nrow = n.traps, ncol = n.traps)
            diag(cov) <- sigma.u
        } else if (cov.structure == "exponential"){
            ## Extracting parameters.
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            rho <- cov.pars$rho
            ## Specifying covariance.
            cov <- sigma.u^2*exp(-trap.dists/rho)
        } else if (cov.structure == "matern"){
            stop("Matern covariance not yet implemented.")
        } else if (cov.structure == "individual"){
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            cov <- matrix(sigma.u^2, nrow = n.traps, ncol = n.traps)    # all covariance elements equal to variance?
        } else if (cov.structure == "lc_exponential"){
            stop("Linear combination of exponentials not yet implemented.")
        } else if (cov.structure == "sq_exponential"){
            ## Extracting parameters.
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            rho <- cov.pars$rho
            ## Specifying covariance.
            cov <- sigma.u^2*exp(-(trap.dists^2)/(rho^2))
        }
        
        ## Simulating random effects.
        u.mat <- mvtnorm::rmvnorm(n, rep(mu.u, n.traps), cov)
        
        ## Getting full encounter rates
        full.er <- base.er*exp(u.mat)
        if (output == 'umat')
            apply(exp(u.mat),1,sum)
        else if (output == 'hazard') 
            apply(full.er,1,sum)
        else if (output == 'pdot')
            1 - exp(-apply(full.er,1,sum))
        else if (output == 'capt') { # capt 
            full.prob <- 1 - exp(-full.er)
            ## Generating Poisson capture histories.
            capt <- matrix(rpois(n*n.traps, full.er), nrow = n)
            capt <- capt[apply(capt, 1, function(x) sum(x) > 0), ]
            convertCH (capt, traps)
        }
        else stop ("unrecognised output type")
    }
}

sim.ou <- function(mu, tau, sigma, n.steps, start = NULL){
    if (is.null(start)){
        start <- mvtnorm::rmvnorm(1, mu, sigma^2*diag(2))
    }
    ## Changing Theo's parameterisation.
    b <- -1/tau
    v <- sigma^2
    out <- matrix(0, nrow = n.steps, ncol = 2)
    out[1, ] <- start
    for (i in 2:n.steps){
        out[i, ] <- mvtnorm::rmvnorm(1, mu + exp(b)*(out[i - 1, ] - mu),
                            v*(1 - exp(2*b))*diag(2))
        
    }
    out
}

## Function to count the number of detections by each trap.
count.dets <- function(locs, traps, epsilon){
    # loc.dists <- crossdist(locs[, 1], locs[, 2],
    #                        traps[, 1], traps[, 2])
    loc.dists <- secr::edist(locs,traps)
    apply(loc.dists, 2, function(x, epsilon) sum(x <= epsilon), epsilon = epsilon)
}

## Closure for detection function.
detfn.closure <- function(detfn, pars){
    if (detfn == "hn"){
        g0 <- pars$g0
        sigma <- pars$sigma
        out <- function(d, er = FALSE){
            out <- g0*exp(-d^2/(2*sigma^2))
            if (er){
                out <- -log(1 - out)
            }
            out
        }
    } else if (detfn == "hr"){
        g0 <- pars$g0
        sigma <- pars$sigma
        z <- pars$z
        out <- function(d, er = FALSE){
            out <- g0*(1 - exp(-((d/sigma)^-z)))
            if (er){
                out <- -log(1 - out)
            }
        }
    } else if (detfn == "hhn"){
        lambda0 <- pars$lambda0
        sigma <- pars$sigma
        out <- function(d, er = FALSE){
            out <- lambda0*exp(-d^2/(2*sigma^2))
            if (!er){
                out <- 1 - exp(-out)
            }
            out
        } 
    } else if (detfn == "hhr"){
        lambda0 <- pars$lambda0
        sigma <- pars$sigma
        z <- pars$z
        out <- function(d, er = FALSE){
            out <- -lambda0*(1 - exp(-((d/sigma)^-z)))
            if (!er){
                out <- 1 - exp(-out)
            }
            out
        }
    }
    out
}

## Setting up detector locations.
trps <- make.grid(4,4,detector = 'count', spacing = 100)
trps <- rbind(trps,trps)
trps$x <- trps$x + rep(c(-10,10), c(16,16))
msk <- make.mask(trps, buffer = 500, type='trapbuffer')

trps2 <- make.grid(6,6,detector='count', spacing = 60)
msk2 <- make.mask(trps2, buffer = 500, type='trapbuffer')

plot(trps)
plot(trps2,add=T, detpar=list(pch=16))
 
