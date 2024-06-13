# 2024-06-07
# Moran.R
# shared code for RSF and SARE

# extractfn for raw capthist 
# autocorrelation of X covariate at detectors, 
# assumed to be saved as attribute 'X' of ch.
# Weight window considers adjacent detectors (queen's move).
MoranX <- function (ch, varname = 'X', ...) {
    if (is.null(attr(ch, varname))) stop (var, " attribute not found")
    X <- attr(ch, varname)
    Xraster <- raster(as.mask(traps(ch)), values = X)
    out <- c(
        n          = nrow(ch),
        detections = sum(ch),
        dnpa       = sum(ch)/nrow(ch), 
        rpsv       = RPSV(ch, CC = TRUE),
        mean.X     = mean(X),
        median.X   = median(X),
        sd.X       = sd(X),
        MoranI     = raster::Moran(Xraster, ...)
    )
    names(out)[5:7] <- paste0(c('mean.', 'median.', 'sd.'), varname)
    out
}

MoranXSummary <- function(sims, var = 'rpsv') {
    output <- summary(sims)$OUTPUT
    var <- match(var, rownames(output[[1]]))
    sumX <- sapply(output, '[', var, c('mean','se'))
    apply(sumX,1,unlist)
}
