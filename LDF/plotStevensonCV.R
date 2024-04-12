library(secr)

# load simulation code and trps, msk
source('D:/Density secr 4.6/knitr/secr-simulations/Stevensontrial.r')

one <- function (r,
                    popn,
                    traps, 
                    det.pars = list(lambda0 = 0.5, sigma = 100),
                    cov.pars = list(sigma.u = 1.2, rho = 150),
                    cov.structure = "sq_exponential",
                    ...) {
    sim.sscr.mge (
        popn     = popn,
        traps    = traps, 
        cov.structure = cov.structure,
        det.pars = det.pars, 
        cov.pars = cov.pars,
        ...)
}

pd <- pdot(msk, trps, detectfn='HHN', detectpar=list(lambda0 = 0.5, sigma = 100), noccasions = 1)
pd2 <- pdot(msk2, trps2, detectfn='HHN', detectpar=list(lambda0 = 0.5, sigma = 100), noccasions = 1)

hd <- -log(1-pd)
one.u.CV <- function (sigma.u = 1.2, rho = 150, mask = msk, traps = trps,
                      covstr = "sq_exponential", nrepl = 500, 
                      output = c('pdot', 'hazard'), ex = 1, ...) {
    output <- match.arg(output)
    tmp <- sapply(1:nrepl, one, 
                  mask, 
                  traps, 
                  cov.pars = list(sigma.u = sigma.u, rho = rho),
                  output=output)
    cv <- apply(tmp, 1, CV)
    pd <- pdot(mask, traps, detectfn='HHN', 
               detectpar=list(lambda0 = 0.5, sigma = 100), noccasions = 1)
    if (output == 'hazard') pd <- -log(1-pd)
        sum(cv * pd^ex) / sum(pd^ex)
}
cvwt <- sapply(seq(0,1.6,0.2), one.u.CV)        
cvwt25 <- sapply(seq(0,1.6,0.2), one.u.CV, ex=0.25)
cvwt25.0 <- sapply(seq(0,1.6,0.2), one.u.CV, ex=0.25, rho=1)
cvwt25.i <- sapply(seq(0,1.6,0.2), one.u.CV, ex=0.25,  covstr="individual")

out <- readRDS(file = 'simldf.RDS')
out0rho <- readRDS(file = 'simldf0rho.RDS')
outi <- readRDS(file = 'simldfi.RDS')

plot(cvwt25, sapply(out, apply, 2, mean, na.rm=TRUE)[1,])
points(cvwt25.0, sapply(out0rho, apply, 2, mean, na.rm=TRUE)[1,], pch=16)
points(cvwt25.i, sapply(outi, apply, 2, mean, na.rm=TRUE)[1,], pch=3)

tmp <- sapply(1:1000, one, msk, trps, output='pdot', cov.structure = 'individual')
tmpind <- sapply(1:1000, one, msk, trps, output='pdot', cov.structure = 'independent')
tmpind0 <- sapply(1:1000, one, msk, trps, output='pdot', cov.pars = list(sigma.u = 1.2, rho = 1))

tmpind0.45 <- sapply(1:1000, one, msk, rotate(trps,45), output='pdot', cov.structure = 'independent')

tmpumat <- as.numeric(sapply(1, one, msk, trps, output='umat'))
hist(log(tmpumat))

covariates(msk) <- data.frame(
    mn    = apply(tmp,1,mean), 
    cv    = apply(tmp,1,CV), 
    cvtrunc = pmin(2, apply(tmp,1,CV)), 
    cvind = apply(tmpind,1,CV),
    cvindtrunc = pmin(2, apply(tmpind,1,CV)), 
    cvind0trunc = pmin(2, apply(tmpind0,1,CV)), 
    cvind045trunc = pmin(2, apply(tmpind0.45,1,CV)), 
    pd    = pd, 
    cvwtd = apply(tmp,1,CV) * pd, 
    umat  = tmpumat)

tmp2 <- sapply(1:1000, one, msk2, trps2, output='pdot')
covariates(msk2) <- data.frame(
    cv2    = apply(tmp2,1,CV), 
    cv2trunc = pmin(2, apply(tmp2,1,CV)),
    pd2    = pd2, 
    cv2wtd = apply(tmp2,1,CV) * pd2 
)

par(mfrow=c(1,3))
plot(msk, cov='pd', dots=F)
plot(msk, cov='mn', dots=F)
plot(msk, cov='cvtrunc', dots=F); plot(trps, add=T)
plot(msk, cov='cvindtrunc', dots=F); plot(trps, add=T)
plot(msk, cov='cvind045trunc', dots=F); plot(rotate(trps,45), add=T)

plot(msk, cov='pd', dots=F)
plot(msk, cov='cvwtd', dots=F)
plot(msk, cov='umat', dots=F)
plot(msk, cov='SD', dots=F)

sum(covariates(msk)$cvwtd) / sum(pd)
sum(covariates(msk2)$cv2wtd) / sum(pd2)
# [1] 0.2769937

par(mfrow=c(1,3))

simcv <- function(sigma.u = 1.2, plt = TRUE, output='pdot') {
    tmp <- sapply(1:200, one, msk, trps, 
                  cov.pars = list(sigma.u = mu, rho = 150), output=output)
    covariates(msk) <- data.frame(cv = apply(tmp,1,CV), pd = pd, cvwtd = apply(tmp,1,CV) * pd)
    if (plt) {
        plot(msk, cov='cv', dots=F)
        plot(msk, cov='pd', dots=F)
        plot(msk, cov='cvwtd', dots=F)
    }
    sum(covariates(msk)$cvwtd) / sum(pd)
}

par(mfrow=c(1,1))
cv <- sapply(seq(0,1.6,0.2), simcv, plt = FALSE, hazard = TRUE)
plot(seq(0,1.6,0.2), cv)
abline(v=1.2, lty=2)

par(mfrow=c(1,3))
cv <- sapply(1.2, simcv, plt = TRUE, hazard = TRUE)

