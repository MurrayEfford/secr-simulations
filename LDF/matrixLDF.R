library(secrdesign)
setNumThreads(24)
nrepl <- 200

simLDF.capthist <- function (traps, popn, detectfn, detectpar, noccasions = 1, 
                             cov.structure, cov.pars, shared = FALSE, ...)
{
    K <- nrow(traps)
    N <- nrow(popn)
    ac.dists   <- secr::edist(popn, traps)     # activity centres to traps
    trap.dists <- secr::edist(traps,traps)     # between traps
    
    # "baseline" encounter rates
    if (detectfn == 14 || detectfn == "HHN"){
        base.er <- detectpar$lambda0*exp(-ac.dists^2/(2*detectpar$sigma^2))
    } 
    else if (detectfn == 16 || detectfn == "HEX"){
        base.er <- detectpar$lambda0*exp(-ac.dists/detectpar$sigma)
    } 
    else stop ("detectfn not implemented")
    
    ## Constructing covariance matrix
    sigma.u <- cov.pars$sigma.u
    rho <- cov.pars$rho
    # always adjust mu.u so E(exp(u.mat[i,j])) = 1.0
    mu.u <- - sigma.u^2 / 2
    if (cov.structure == "none"){
        cov <- matrix(0, nrow = K, ncol = K)
    } 
    else if (cov.structure == "independent"){    # rho = 0
        cov <- diag(sigma.u^2, K, K)             
    } 
    else if (cov.structure == "exponential"){
        cov <- sigma.u^2*exp(-trap.dists/rho)
    } 
    else if (cov.structure == "individual"){     # rho > infty
        cov <- matrix(sigma.u^2, nrow = K, ncol = K)
    } 
    else if (cov.structure == "sq_exponential"){
        cov <- sigma.u^2*exp(-(trap.dists^2)/(rho^2))
    }
    ## Simulate random field (individual x detector)
    if (shared) {
        u.row <- mvtnorm::rmvnorm(1, mean = rep(mu.u, K), sigma = cov) 
        u.mat <- matrix(rep(u.row,N), byrow = TRUE, nrow = N)
    }
    else {    
        u.mat <- mvtnorm::rmvnorm(N, mean = rep(mu.u, K), sigma = cov) 
    }
    full.er <- base.er * exp(u.mat) # encounter rates
    if (noccasions>1) {
        # replicate occasions
        full.er <- aperm(apply(full.er,1:2,rep,noccasions), c(2,1,3))
    }
    if (detector(traps)[1] == 'proximity') {
        ch <- rbinom(N*K*noccasions, 1, 1-exp(-full.er))
    }
    else if (detector(traps)[1] == 'count') {
        ch <- rpois(N*K*noccasions, full.er)
    }
    else stop ("detector not proximity or count")
    ch <- array(ch, dim=c(N,noccasions,K))
    rownames(ch) <- 1:N
    ch <- ch[apply(ch, 1, function(x) sum(x) > 0),,, drop = FALSE]
    class(ch) <- 'capthist'
    traps(ch) <- traps
    attr(ch, 'u.mat') <- u.mat
    ch
}
trps <- make.grid(4,4,detector = 'count', spacing = 100)
trps <- rbind(trps,trps)
trps$x <- trps$x + rep(c(-10,10), c(16,16))
msk    <- make.mask(trps, buffer = 500, type = 'trapbuffer')
expand.detarg <- function (sigma.u, rho, ...) {
    arglist <- list(...)
    comb <- expand.grid(sigma.u = sigma.u, rho = rho, KEEP.OUT.ATTRS = FALSE)
    arglist <- rep(list(arglist), nrow(comb))
    for (i in 1:length(arglist)) arglist[[i]]$cov.pars <- as.list(comb[i,])
    arglist
}
detarg <- expand.detarg(sigma.u = seq(0,1.6,0.4), rho = c(50,150,450), 
                        cov.structure = "sq_exponential", noccasions = 1)
poparg <- list(buffer = 500)
fitarg <- list(detectfn = 'HHN', mask = msk)
scen   <- make.scenarios(D = 1, lambda0 = 0.5, sigma = 100,  detectfn = 'HHN', 
                         detindex = 1:length(detarg))

LDFsimsold2 <- run.scenarios(nrepl, scen, trapset = trps, CH.function = "simLDF.capthist",
                             pop.args = poparg, det.args = detarg, fit = TRUE, fit.args = fitarg)
saveRDS(LDFsimsold2, file = 'LDFsimsold2.RDS')

detarg <- expand.detarg(sigma.u = seq(0,1.6,0.4), rho = c(50,150,450), 
                        cov.structure = "sq_exponential", shared = TRUE, noccasions = 1)
poparg <- list(buffer = 500)
fitarg <- list(detectfn = 'HHN', mask = msk)
scen   <- make.scenarios(D = 1, lambda0 = 0.5, sigma = 100,  detectfn = 'HHN', 
                         detindex = 1:length(detarg))

LDFsimsold2shared <- run.scenarios(nrepl, scen, trapset = trps, CH.function = "simLDF.capthist",
                                   pop.args = poparg, det.args = detarg, fit = TRUE, fit.args = fitarg)
saveRDS(LDFsimsold2shared, file = 'LDFsimsold2shared.RDS')
