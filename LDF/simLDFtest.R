library(secr)

grid <- make.grid(spacing=40, detector = 'proximity')

one <- function(r) {
    pop <- sim.popn(D=10, grid)
    ch <- simLDF.capthist (grid, pop, 14, detectpar = list(lambda0=0.5, sigma=20), noccasions = 5,
                           cov.structure = "exponential", cov.pars = list(sigma.u=1, rho = 200), 
                           shared = FALSE)
    message("mean(exp(u.mat)) ", mean(exp(attr(ch,'u.mat'))))
    predict(secr.fit(ch, detectfn = 'HHN', buffer = 100, trace = FALSE))
    # mean(exp(attr(ch,'u.mat')))
}

out<- lapply(1:100, one)

mean(sapply(out,'[','D','estimate')) / 10 -1
mean(sapply(out,'[','lambda0','estimate')) / 0.5 -1
mean(sapply(out,'[','sigma','estimate')) / 0.5 -1

