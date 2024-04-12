## ORNSTEIN-UHLENBECK SIMULATION CF STEVENSON ET AL 2021 USING sscr CODE
## 2024-03-12 etc.

library(secr)
source('D:/Density secr 4.6/knitr/secr-simulations/Stevensontrial.r')

# Demo OU from Stevenson code (modified from Theo)
sig <- 100
par(mfrow=c(2,2), pty='s')
MASS::eqscplot(sim.ou(mu = c(0,0), tau = 0.1, sigma = sig, n.steps=1000), type='o')
MASS::eqscplot(sim.ou(mu = c(0,0), tau = 1,   sigma = sig, n.steps=1000), type='o')
MASS::eqscplot(sim.ou(mu = c(0,0), tau = 10,  sigma = sig, n.steps=1000), type='o')
MASS::eqscplot(sim.ou(mu = c(0,0), tau = 100, sigma = sig, n.steps=1000), type='o')

checkfn <- function(tau = 1, sigma = 100, nsteps = 10000) {
    sims <- sim.ou(mu = c(0,0), tau = tau, sigma = sigma, n.steps=nsteps)
    c(apply(sims,2,mean), apply(sims,2,sd) )
}
checkfn()


## Simulating data.
simfnOU <- function (nrepl = 1,
                    traps, 
                    mask, 
                    D        = 1, 
                    det.pars = list(lambda0 = 0.5, sigma = 100),
                    cov.pars = list(),
                    resp     = "pois",
                    cov.structure = "sq_exponential",
                    ...) {
    out <- vector(mode = 'list', nrepl)
    for (r in 1:nrepl) {
        # simplified function from Stevensontrial.r
        pop <- sim.popn(D, traps, buffer = 500)
        ch <- sim.sscr.mge(
            popn     = pop, 
            traps    = traps, 
            cov.structure = cov.structure,
            det.pars = det.pars, 
            cov.pars = cov.pars,
            output   = 'capt')
        pred <- predict(secr.fit(ch, mask = msk, ...))
        out[[r]] <- list(pred=pred, n = nrow(ch), detections=sum(ch))
    }
    out
}

setNumThreads(18)
set.seed <- 123
one.OU <- function (tau = 10, 
                   cov.pars = list(sigma=100, n.steps = 100, epsilon = 10),
                   nrepl = 500,
                   ...) {
    cov.pars$tau <- tau
    simfnOU(nrepl, traps = trps, mask = msk,
           cov.structure = "OU", cov.pars = cov.pars,
           detectfn = 'HHN', trace = FALSE, ...)
}

## test
sims <- one.OU(tau=100, nrepl=10)
getD <- function (sims, parm = 'D') mean(sapply(sims, function(x) x$pred[parm,'estimate']))
getD(sims)

## run
OUout <- lapply(c(1,10,100), one.OU, nrepl = 100)
sapply(OUout, getD)

saveRDS(OUout, file = 'simOU.RDS')

sapply(OUout, apply, 2, mean, na.rm=T)

x <- seq(0,1.6,0.2)
OUout <- readRDS(file = 'simOU.RDS')

par(mfrow=c(1,3), mar=c(5,5,2,2))
plot(x, sapply(OUout, apply, 2, mean, na.rm=TRUE)[1,], ylim=c(0.4,1.2), xlab='sigma.u',ylab='D-hat')
abline(h=1)


## capthist simulation function that can be passed to secrdesign::run.scenarios
simOU.capthist <- function (
        traps,
        popn,
        detectpar = list(),
        noccasions = 100,   # effective "duration"
        epsilon = 10,    # proximity threshold for detection
        ...)
{
    # adapted from Supplements of Ben Stevenson 2021 Biometrics paper
    sim.ou <- function(mu, tau, sigma, n.steps){
        start <- mvtnorm::rmvnorm(1, mu, sigma^2*diag(2))
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
    count.dets <- function(locs, traps, epsilon){
        loc.dists <- secr::edist(locs,traps)
        apply(loc.dists, 2, function(x, epsilon) sum(x <= epsilon), epsilon = epsilon)
    }

    # re-purpose parameters:
    #    detectpar$lambda0 as tau
    #    noccasions as n.steps
    n <- nrow(popn)
    tau <- detectpar$lambda0
    sigma <- detectpar$sigma
    n.steps <- noccasions
    capt <- matrix(0, nrow = n, ncol = nrow(traps))
    for (i in 1:n){
        mu <- as.numeric(popn[i, ])   # AC
        locs <- sim.ou(mu, tau, sigma, n.steps)
        capt[i, ] <- count.dets(locs, traps, epsilon)
    }
    capt <- capt[apply(capt, 1, function(x) sum(x) > 0), ]
    ch <- array(capt, dim=c(n,1,ncol(capt)), dimnames=list(1:n,NULL,NULL))
    class(ch) <- 'capthist'
    traps(ch) <- traps
    ch
}


library(secrdesign)
scen <- make.scenarios(D = 1, lambda0 = c(1,10,100), sigma = 100, noccasions = 100, detectfn = 'HHN')
testdet <- run.scenarios(scen, nrepl = 10, trapset = trps, 
                         CH.function = "simOU.capthist",
                         det.args = list(epsilon = 10))
summary(testdet)

testfit <- run.scenarios(scen, 
                         nrepl       = 100, 
                         trapset     = trps, 
                         fit         = TRUE,
                         CH.function = "simOU.capthist",
                         det.args    = list(epsilon = 10),
                         fit.args    = list(detectfn = 'HHN'))
summary(testfit)