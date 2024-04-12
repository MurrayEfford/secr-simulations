# function to return simulated capthist object for 
# variable, autocorrelated detector efficiency VAC

simVAC.capthist <- function (traps, popn, detectfn, detectpar, noccasions = 1, 
                             betaX = -0.5, phi = 1, cov.structure = "exponential",
                             usermvn = FALSE, scale = FALSE, makeuniform = FALSE,...)
{
    rmvn <- function(n, mu = 0, V = matrix(1), seed = NULL) {
        # based on Moqanaki 2021 code
        p <- length(mu)
        D <- V   # no Cholesky decomposition D <- chol(V)
        t(matrix(rnorm(n * p, mean = 0, sd = 1), ncol = p) %*% D + rep(mu, rep(n, p)))
    }
    
    K <- nrow(traps)
    N <- nrow(popn)
    trap.dists <- secr::edist(traps,traps)     # between traps
    
    if (cov.structure == "exponential")
        cov <- exp(-phi * trap.dists)     
    else if (cov.structure == "sq_exponential")
        cov <- exp(-phi^2 * trap.dists^2)     
        
    if (usermvn) 
        X.row <- rmvn(1, mu = rep(0, K), V = cov)   # as in Moqanaki et al.
    else 
        X.row <- mvtnorm::rmvnorm(1, mean = rep(0, K), sigma = cov) 
    
    X.row <- as.numeric(X.row)
    if (scale) X.row <- scale(X.row)

    if (makeuniform) {
        # this code transforms the MVN values so they are uniformly distributed on [-1.96,1.96].
        df <- data.frame(row = 1:K, orig.value = X.row)
        df <- df[order(df$orig.value),]
        df$value <- seq(-1.96, 1.96, length.out = K)
        df <- df[order(df$row), ]
        X.row <- df$value
    }

    # replicate across animals (rows)
    X.mat <- matrix(X.row, byrow = TRUE, nrow = N, ncol = K)
    
    # Moqanaki et al. (2021) Eq 8,9
    beta0 <- logit(detectpar$g0)   # intercept when no variation
    p0 <- invlogit(beta0 + betaX * X.mat)
    
    # distance effect on detection probability
    ac.dists   <- secr::edist(popn, traps)     # activity centres to traps
    if (detectfn == 0 || detectfn == "HN"){
        p <- p0 * exp(-ac.dists^2/(2*detectpar$sigma^2))
    } 
    else stop ("detectfn not implemented")
                   
    if (noccasions>1) {
        p <- aperm(apply(p,1:2,rep,noccasions), c(2,1,3))
    }
    if (detector(traps)[1] == 'proximity') {
        ch <- rbinom(N*K*noccasions, 1, p)
    }
    else if (detector(traps)[1] == 'count') {
        ch <- rpois(N*K*noccasions, -log(1-p))
    }
    else stop ("detector not proximity or count")
    ch <- array(ch, dim=c(N,noccasions,K))
    rownames(ch) <- 1:N
    ch <- ch[apply(ch, 1, function(x) sum(x) > 0),,, drop = FALSE]
    class(ch) <- 'capthist'
    traps(ch) <- traps
    attr(ch, 'X.row') <- as.numeric(X.row)
    attr(ch, 'p0') <- invlogit(beta0 + betaX * X.row)
    ch
}

grid400 <- make.grid(20, 20, spacex = 1, detector = 'proximity')
mask400 <- make.mask(grid400, buffer = 4.5)
pop <- sim.popn(Nbuffer = 250, core = grid400, Ndist = 'fixed', buffer = 4.5)

ch <- simVAC.capthist (grid400, pop, detectfn=0, detectpar=list(g0 = 0.15, sigma = 1.5),
                       noccasions = 1, betaX = -0.5, phi = 0.001, usermvn = FALSE, scale=F, 
                       makeuniform = TRUE, cov.structure="sq_exponential")
summary(ch)

XMoran <- function (ch) {
    Xraster <- raster(as.mask(traps(ch)), values = attr(ch, 'X.row'))
    raster::Moran(Xraster)
}
XMoran(ch)

mean(attr(ch, 'p0'))
CV(attr(ch, 'p0'))

covariates(grid400) <- data.frame(X = attr(ch, 'X.row'))
plot(as.mask(grid400), cov = 'X', dots = FALSE)

#################################

# test alg

testfn <-  function (r, traps = grid400, phi = 1, usermvn = FALSE, scale = FALSE, 
                     makeuniform = FALSE, plt = FALSE, raw = FALSE, sqexp = FALSE...) {
    rmvn <- function(n, mu = 0, V = matrix(1), seed = NULL) {
        p <- length(mu)
        if (any(is.na(match(dim(V), p))))
            stop("Dimension problem!")
        #D <- chol(V)
        D <- V
        if(!is.null(seed)){set.seed(seed)}
        t(matrix(rnorm(n * p, mean = 0, sd = 1), ncol = p) %*% D + rep(mu, rep(n, p)))
    }
    trap.dists <- secr::edist(traps,traps)     # between traps
    K <- nrow(traps)
    cov <- exp(-phi * trap.dists)
    if (sqexp) cov <- cov%*%t(cov)
    
    if (usermvn) {
        X <- rmvn(1, mu = rep(0, K), V = cov) 
    }
    else {
        X <- mvtnorm::rmvnorm(1, mean = rep(0, K), sigma = cov) 
    }
    X <- as.numeric(X)
    if (scale) X <- scale(X)
    if (makeuniform) {
        df <- data.frame(row = 1:K, orig.value = X)
        df <- df[order(df$orig.value),]
        df$value <- seq(-1.96, 1.96, length.out = K)
        df <- df[order(df$row), ]
        X <- df$value
    }
    if (raw) {
        return (list(cov = as.numeric(cov), X=X, dists = as.numeric(trap.dists)))
    }
    else {
        Xraster <- raster::rasterFromXYZ(cbind(traps, X))
        if (plt) raster::plot(Xraster)
        raster::Moran(Xraster)
    }
}
testfn(1, phi = 0.001, sqexp = TRUE)

t001 <- sapply(1:1000, testfn, phi=0.001)
t001m <- sapply(1:1000, testfn, phi=0.001, usermvn = FALSE, makeuniform = TRUE)
t001u <- sapply(1:1000, testfn, phi=0.001, usermvn = TRUE, makeuniform = FALSE)
t001um <- sapply(1:1000, testfn, phi=0.001, usermvn = TRUE, makeuniform = TRUE)
t001ums <- sapply(1:1000, testfn, phi=0.001, usermvn = TRUE, makeuniform = TRUE, scale = TRUE)
t001q <- sapply(1:1000, testfn, phi=0.001, sqexp = TRUE)

sapply(list(t001, t001m, t001u, t001um, t001ums), median)
# [1] 0.83041 0.82725 0.94030 0.95386 0.95321
sapply(list(t001, t001m, t001u, t001um, t001ums), sd)/1000^0.5
# [1] 0.00219931 0.00229788 0.00070886 0.00046434 0.00044562

#################################

library(secrdesign)
poparg <- list(Ndist = 'fixed')
detarg <- expand.arg(phi = c(0.001,1,1000), betaX = c(-0.5,-2), usermvn = c(FALSE,TRUE))
parm <- do.call(rbind, lapply(detarg, as.data.frame))   # now flat
parm

fitarg <- list(details = list(distribution = 'binomial'))
scen <- make.scenarios(noccasions = 1, D=250/maskarea(mask400), g0 = 0.15, 
                       sigma = 1.5, detectfn = 0, detindex = 1:length(detarg))
checkMoran <- run.scenarios(10, scen, grid400, maskset=list(mask400), CH.function = "simVAC.capthist", 
         pop.args = poparg, det.args=detarg, fit = FALSE, extractfn = XMoran)
# not makeuniform
#      phi betaX usermvn     MoranI       min      max
# 1  1e-03  -0.5   FALSE  0.8062388  0.722335 0.942718
# 2  1e+00  -0.5   FALSE  0.2910467  0.213978 0.374047
# 3  1e+03  -0.5   FALSE  0.0155227 -0.022388 0.051939
# 4  1e-03  -2.0   FALSE  0.8852857  0.598274 0.904026
# 5  1e+00  -2.0   FALSE  0.2853841  0.205200 0.309699
# 6  1e+03  -2.0   FALSE  0.0020345 -0.031853 0.020615
# 7  1e-03  -0.5    TRUE  0.9432546  0.896146 0.956395
# 8  1e+00  -0.5    TRUE  0.6545223  0.555853 0.718077
# 9  1e+03  -0.5    TRUE  0.0027837 -0.041145 0.069033
# 10 1e-03  -2.0    TRUE  0.9246241  0.895132 0.957289
# 11 1e+00  -2.0    TRUE  0.6162806  0.558712 0.722437
# 12 1e+03  -2.0    TRUE -0.0051913 -0.066712 0.060017

data.frame(parm, 
           MoranI =sapply(checkMoran$output, median), 
           min =sapply(checkMoran$output, min),
           max =sapply(checkMoran$output, max))


detarg <- expand.arg(phi = c(0.001,1,1000), betaX = c(-0.5,-2), usermvn = TRUE, makeuniform = TRUE)
fitarg <- list(details = list(distribution = 'binomial'))
scen <- make.scenarios(noccasions = 1, D=250/maskarea(mask400), g0 = 0.15, 
                       sigma = 1.5, detectfn = 0, detindex = 1:length(detarg))

simsVAC <- run.scenarios(100, scen, grid400, maskset=list(mask400), CH.function = "simVAC.capthist", 
                         pop.args = poparg, det.args=detarg, fit = TRUE, fit.args = fitarg)

estimateSummary(simsVAC)
#   scenario true.D nvalid    EST  seEST         RB      seRB      RSE    RMSE    rRMSE  COV
# 1        1 3188.8    100 2914.4 22.587 -0.0860565 0.0070833 0.064235  354.70 0.111233 0.64
# 2        2 3188.8    100 3109.8 22.200 -0.0247816 0.0069620 0.064761  234.60 0.073570 0.90
# 3        3 3188.8    100 3196.4 21.583  0.0023842 0.0067685 0.064863  214.89 0.067388 0.95
# 4        4 3188.8    100 2173.8 23.964 -0.3183070 0.0075151 0.052901 1042.64 0.326972 0.01
# 5        5 3188.8    100 2882.2 21.208 -0.0961466 0.0066510 0.049794  372.19 0.116720 0.51
# 6        6 3188.8    100 3202.3 14.119  0.0042448 0.0044278 0.047993  141.13 0.044260 0.97

detarg1 <- expand.arg(phi = c(0.001,1,1000), betaX = c(-0.5,-2), usermvn = FALSE, makeuniform = TRUE)
simsVAC1 <- run.scenarios(100, scen, grid400, maskset=list(mask400), CH.function = "simVAC.capthist", 
                         pop.args = poparg, det.args=detarg1, fit = TRUE, fit.args = fitarg)
do.call(rbind, lapply(detarg1, as.data.frame))
#     phi betaX usermvn makeuniform
# 1 1e-03  -0.5   FALSE        TRUE
# 2 1e+00  -0.5   FALSE        TRUE
# 3 1e+03  -0.5   FALSE        TRUE
# 4 1e-03  -2.0   FALSE        TRUE
# 5 1e+00  -2.0   FALSE        TRUE
# 6 1e+03  -2.0   FALSE        TRUE

saveRDS(simsVAC, file = 'simsVAC.RDS')    # usermvn
saveRDS(simsVAC1, file = 'simsVAC1.RDS')  # rmvtnorm

estimateSummary(simsVAC1)
#   scenario true.D nvalid    EST  seEST         RB      seRB      RSE   RMSE    rRMSE  COV
# 1        1 3188.8    100 2989.6 23.273 -0.0624712 0.0072984 0.064346 305.46 0.095792 0.73
# 2        2 3188.8    100 3158.1 20.505 -0.0096223 0.0064304 0.064896 206.32 0.064701 0.95
# 3        3 3188.8    100 3196.4 21.583  0.0023842 0.0067685 0.064863 214.89 0.067388 0.95
# 4        4 3188.8    100 2439.7 21.032 -0.2348998 0.0065956 0.051536 777.73 0.243895 0.02
# 5        5 3188.8    100 3083.2 17.609 -0.0331109 0.0055222 0.048712 204.56 0.064150 0.87
# 6        6 3188.8    100 3202.3 14.119  0.0042448 0.0044278 0.047993 141.13 0.044260 0.97

estimateSummary(simsVAC1, 'g0')
#   scenario true.g0 nvalid     EST     seEST       RB     seRB      RSE     RMSE   rRMSE  COV
# 1        1    0.15    100 0.17998 0.0017432 0.199843 0.011622 0.100062 0.034633 0.23089 0.57
# 2        2    0.15    100 0.16844 0.0016273 0.122915 0.010849 0.102058 0.024538 0.16359 0.83
# 3        3    0.15    100 0.16375 0.0017861 0.091669 0.011907 0.102314 0.022470 0.14980 0.86
# 4        4    0.15    100 0.42424 0.0039680 1.828276 0.026454 0.059000 0.277069 1.84713 0.00
# 5        5    0.15    100 0.32010 0.0025616 1.133977 0.017077 0.063351 0.171996 1.14664 0.00
# 6        6    0.15    100 0.29606 0.0016475 0.973735 0.010983 0.063249 0.146977 0.97985 0.00

estimateSummary(simsVAC1, 'sigma')
#   scenario true.sigma nvalid    EST     seEST         RB      seRB      RSE     RMSE    rRMSE  COV
# 1        1        1.5    100 1.4815 0.0055460 -0.0123369 0.0036973 0.041485 0.058202 0.038801 0.95
# 2        2        1.5    100 1.4878 0.0059744 -0.0081609 0.0039829 0.042133 0.060691 0.040461 0.94
# 3        3        1.5    100 1.5046 0.0061473  0.0030482 0.0040982 0.042276 0.061336 0.040891 0.96
# 4        4        1.5    100 1.4224 0.0039716 -0.0517163 0.0026477 0.025888 0.087059 0.058040 0.48
# 5        5        1.5    100 1.4629 0.0048999 -0.0247360 0.0032666 0.027083 0.061267 0.040845 0.77
# 6        6        1.5    100 1.4977 0.0038673 -0.0015267 0.0025782 0.026989 0.038547 0.025698 0.95


detarg2 <- expand.arg(phi = c(0.001,1,1000), betaX = c(-0.5,-2), usermvn = FALSE, 
                      cov.structure = "sq_exponential", makeuniform = TRUE)
simsVAC2 <- run.scenarios(100, scen, grid400, maskset=list(mask400), CH.function = "simVAC.capthist", 
                          pop.args = poparg, det.args=detarg2, fit = TRUE, fit.args = fitarg)

> estimateSummary(simsVAC2)
scenario true.D nvalid    EST  seEST         RB      seRB      RSE    RMSE    rRMSE  COV
1        1 3188.8    100 2920.0 19.777 -0.0842861 0.0062020 0.064261  333.10 0.104461 0.66
2        2 3188.8    100 3166.7 21.595 -0.0069369 0.0067723 0.064863  216.01 0.067740 0.94
3        3 3188.8    100 3196.4 21.583  0.0023842 0.0067685 0.064863  214.89 0.067388 0.95
4        4 3188.8    100 2169.6 13.133 -0.3196277 0.0041184 0.052661 1027.56 0.322244 0.00
5        5 3188.8    100 3116.2 16.726 -0.0227451 0.0052452 0.048545  181.54 0.056930 0.93
6        6 3188.8    100 3202.3 14.119  0.0042448 0.0044278 0.047993  141.13 0.044260 0.97

print(round(estimateSummary(simsVAC2)[,c(1,3,6,8,11)],3), row.names=F)
scenario nvalid     RB   RSE  COV
       1    100 -0.084 0.064 0.66
       2    100 -0.007 0.065 0.94
       3    100  0.002 0.065 0.95
       4    100 -0.320 0.053 0.00
       5    100 -0.023 0.049 0.93
       6    100  0.004 0.048 0.97
