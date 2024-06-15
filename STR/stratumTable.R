# stratumTable.R

# report subregion density estimates and overall

stratumTable <- function (fit, mask = NULL, cov = 'stratum', poly = NULL, alpha = 0.05) {
    if (is.null(mask)) mask <- fit$mask
    if (!is.null(poly)) mask <- subset(mask, pointsInPolygon(mask, poly))
    if (is.null(covariates(mask)) || !cov %in% names(covariates(mask))) {
        warning("covariate missing from mask, assuming equal areas")
        if (inherits(fit,'secr')) {
            pred <- predict(fit, all.levels = TRUE)
        }
        else {
            pred <- fit   # assume predicted values
        }
        Ah <- rep(1,length(pred))
        strata <- factor(1:length(pred))
    }
    else {
        strata <- covariates(mask)[,cov]
        if (!is.factor(strata)) strata <- factor(strata)
        # if (length(levels(strata))<2) stop ("requires multiple strata")
        freq <- table(strata)
        Ah <- as.numeric(freq * attr(mask, 'area'))
        newdat <- data.frame(stratum = factor(levels(strata)))
        names(newdat)[names(newdat)=='stratum'] <- cov
        pred <- predict(fit, newdata = newdat)
    }
    Wh <- Ah/sum(Ah)
    tab <- data.frame (
        Ah = Ah,
        Wh = Wh,
        D = sapply(pred, '[', 'D','estimate'),
        SED = sapply(pred, '[', 'D','SE.estimate'),
        lcl =  sapply(pred, '[', 'D','lcl'),
        ucl =  sapply(pred, '[', 'D','ucl'),
        row.names = levels(strata)
    )
    Dtot <- sum(tab$D*tab$Wh)
    SEtot <- sqrt( sum(tab$SED^2 * tab$Wh^2 ))
    z <- abs(qnorm(1 - alpha/2))
    total <- data.frame (
        Ah = sum(Ah), 
        Wh = sum(Wh),
        D  = Dtot, 
        SED = SEtot,
        lcl = Dtot - z * SEtot,
        ucl = Dtot + z * SEtot,
        row.names = 'Total')
    rbind(tab, total)
}

# as extractfn...
# pass cov = 'grid'

stratumSummary <- function(sims, trueD = c(5,20,12.5), ...) {
    onescenario <- function (tables) {
        abc <- abind::abind(tables, along=3)
        mn <- as.data.frame(apply(abc, 1:2,mean))
        mn$RB <- (mn$D-trueD)/trueD
        mn$RSE <- mn$SED/mn$D
        mn$COV <- numeric(3)
        for (i in 1:3) mn$COV[i] <- mean (abc[i,'lcl',]<trueD[i] & abc[i,'ucl',]>trueD[i])
        mn
    }
    lapply(sims$output, onescenario)
}
