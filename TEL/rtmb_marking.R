###############################################################################
## rtmb_marking.R
## One-stage maximum likelihood for telemetry-marking data with RTMB, as an independent
## implementation of the model fitted by secr::secr.fit for telemetry type 'marking'.
## Use as the 'extractfn' of secrdesign::run.scenarios(fit = FALSE, ...)
## (see function rtmb_extract at the end of this file).
##
## Model (count detectors, homogeneous density, hazard half-normal, uniform prior for
## the activity centres of the collared animals):
##   * collared animal i: identified cue counts Poisson(q U_k lambda_k(s_i)); n_i telemetry
##     fixes bivariate normal about s_i with the same sigma as the detection function; the
##     activity centre s_i is integrated out by Gauss-Hermite quadrature about the mean fix
##     (or summed over the mask if the animal has no fixes);
##   * unmarked cues: total C_u negative binomial with mean mu = D J1 - sum_i E_i[L(s_i)],
##     where J1 = integral of the cue rate L(x) over the mask and E_i is the posterior mean
##     given animal i's own data, and variance c-hat * mu (compound-Poisson; Poisson or
##     fixed number of animals; the posterior variance of the collared animals' cues is
##     included); c-hat is estimated in a first pass and then held fixed;
##   * unidentified marked cues Tm (optional): Poisson((1 - q) U_k sum_i E_i[lambda_k]).
## Differences from secr: the unmarked total is negative binomial rather than Poisson
## with a variance scale factor (point estimates of D differ by 0.4-2%), and there is no
## density model.
##
## Requires packages RTMB and secr (exported functions only). Units: mask area in
## hectares gives density per hectare.
###############################################################################

## Gauss-Hermite nodes and weights for a standard normal (probabilists'), n nodes
gh_std_normal <- function (n) {
    J <- matrix(0, n, n)
    if (n > 1) for (i in 1:(n - 1)) J[i, i + 1] <- J[i + 1, i] <- sqrt(i)
    e <- eigen(J, symmetric = TRUE)
    list(z = e$values, w = e$vectors[1, ]^2)
}

## Fit to a capthist of telemetry type 'marking'.
## distribution: "binomial" (fixed N) or "poisson". Returns a list, or NULL if the fit failed.
fit_rtmb_marking <- function (capthist, mask, distribution = c("binomial", "poisson"),
                              ngh = 7, trace = FALSE) {
    ## RTMB must be attached (not just loaded): its automatic-differentiation versions of
    ## functions such as plogis and lgamma are used by the likelihood below
    suppressPackageStartupMessages(require(RTMB))
    distribution <- match.arg(distribution)
    ptm <- proc.time()

    ## 1. Data --------------------------------------------------------------
    traps <- secr::traps(capthist)
    ## number of sighting detectors: all rows of 'traps' except the notional telemetry detector
    ## (detector() has one element per occasion, not per detector)
    K0    <- nrow(traps) - any(secr::detector(traps) == "telemetry")
    S     <- sum(secr::markocc(traps) == 0)                   # sighting occasions
    CH0   <- apply(capthist[, 1:S, 1:K0, drop = FALSE], c(1, 3), sum)   # n x K0 identified cues
    traps_xy <- as.matrix(traps[1:K0, ])
    u     <- secr::usage(traps)
    usage <- if (is.null(u)) matrix(1, K0, S) else as.matrix(u)[1:K0, 1:S, drop = FALSE]
    U0    <- rowSums(usage)
    keep  <- U0 > 0                      # unused detectors carry no information
    K     <- sum(keep)
    CH    <- CH0[, keep, drop = FALSE]
    U     <- U0[keep]
    trap_xy <- traps_xy[keep, , drop = FALSE]
    n     <- nrow(CH)

    mask_xy <- as.matrix(mask)
    a       <- attr(mask, "area")
    M       <- nrow(mask_xy)
    A_mask  <- M * a

    Tu  <- attr(capthist, "Tu")
    if (is.null(Tu)) stop("capthist has no Tu attribute")
    C_u <- sum(Tu)
    Tm  <- attr(capthist, "Tm")
    Tm_k <- if (is.null(Tm)) rep(0, K) else rowSums(as.matrix(Tm)[1:K0, 1:S, drop = FALSE])[keep]
    has_Tm <- sum(Tm_k) > 0

    telem <- secr::telemetryxy(capthist)
    has_tel <- if (length(telem) == n) vapply(seq_len(n), function(i)
        !is.null(telem[[i]]) && nrow(telem[[i]]) > 0, logical(1)) else rep(FALSE, n)
    tel_i <- which(has_tel); non_i <- which(!has_tel)
    telinfo <- lapply(tel_i, function(i) {
        tx <- as.matrix(telem[[i]]); tb <- colMeans(tx)
        list(tb = tb, nl = nrow(tx), ssd = sum(sweep(tx, 2, tb)^2)) })

    gh <- gh_std_normal(ngh)
    zz <- as.matrix(expand.grid(z1 = gh$z, z2 = gh$z))
    logw_gh <- log(as.vector(outer(gh$w, gh$w)))
    G  <- nrow(zz)

    d2m <- outer(mask_xy[, 1], trap_xy[, 1], "-")^2 + outer(mask_xy[, 2], trap_xy[, 2], "-")^2  # M x K
    Cm  <- matrix(rep(trap_xy[, 1], each = G), G, K)
    Cy  <- matrix(rep(trap_xy[, 2], each = G), G, K)
    lgY <- rowSums(lgamma(CH + 1))
    dat <- list(CH = CH, U = U, K = K, S = S, n = n, M = M, a = a, A_mask = A_mask,
                d2m = d2m, Cm = Cm, Cy = Cy, G = G, z1 = zz[, 1], z2 = zz[, 2], logw = logw_gh,
                tel_i = tel_i, non_i = non_i, telinfo = telinfo,
                YlogU = as.vector(CH %*% log(U)), sumY = rowSums(CH), lgY = lgY,
                Tm_k = Tm_k, Tm_pos = which(Tm_k > 0), lgTm = sum(lgamma(Tm_k + 1)), has_Tm = has_Tm,
                C_u = C_u, distribution = distribution, chat_on_tape = TRUE, chat = NA_real_)

    ## 2. Starting values ---------------------------------------------------
    s_init <- t(vapply(seq_len(n), function(i) {
        if (has_tel[i]) telinfo[[match(i, tel_i)]]$tb
        else if (any(CH[i, ] > 0)) colMeans(trap_xy[CH[i, ] > 0, , drop = FALSE])
        else colMeans(trap_xy)
    }, numeric(2)))
    sig0 <- if (any(has_tel)) {
        sqrt(sum(vapply(telinfo, function(x) x$ssd, 0)) / (2 * sum(vapply(telinfo, function(x) x$nl, 0))))
    } else secr::RPSV(capthist, CC = TRUE)
    Lam_at <- function(xy, lam0, sig)
        rowSums(sapply(1:K, function(k)
            U[k] * lam0 * exp(-((xy[, 1] - trap_xy[k, 1])^2 + (xy[, 2] - trap_xy[k, 2])^2) / (2 * sig^2))))
    lam00 <- (sum(CH) + sum(Tm_k)) / max(sum(Lam_at(s_init, 1, sig0)), 1e-12)
    lam00 <- min(max(lam00, 1e-6), 10)
    I10   <- sum(Lam_at(mask_xy, lam00, sig0)) * a
    Du0   <- max(C_u, 1) / I10
    q0    <- if (has_Tm) min(max(sum(CH) / (sum(CH) + sum(Tm_k)), 0.1), 0.9) else 1

    pars <- list(log_Du = log(Du0), log_lambda0 = log(lam00), log_sigma = log(sig0),
                 logit_q = qlogis(q0))
    map  <- if (has_Tm) list() else list(logit_q = factor(NA))

    smooth_pos <- function(x, eps) 0.5 * (x + sqrt(x^2 + eps^2))
    lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

    ## 3. Negative log-likelihood ------------------------------------------
    nll_joint <- function (pars) {
        Du      <- exp(pars$log_Du)
        lambda0 <- exp(pars$log_lambda0)
        sigma   <- exp(pars$log_sigma)
        q       <- if (dat$has_Tm) plogis(pars$logit_q) else 1
        two_s2  <- 2 * sigma^2
        K <- dat$K; U <- dat$U

        loglam_m <- pars$log_lambda0 - dat$d2m / two_s2
        lam_m    <- exp(loglam_m)
        Lm  <- (lam_m %*% U)[, 1]                  # cue rate (all occasions) at each mask cell
        J1  <- sum(Lm) * dat$a
        J2  <- sum(Lm^2) * dat$a

        nll <- 0
        SumEL <- 0; Vs <- 0
        Elam  <- rep(0, K)                          # sum_i E_i[lambda_k(s_i)]

        ## collared animals without telemetry: sum over mask cells
        nu <- length(dat$non_i)
        if (nu > 0) {
            Y  <- dat$CH[dat$non_i, , drop = FALSE]
            cu <- dat$YlogU[dat$non_i] + dat$sumY[dat$non_i] * log(q) - dat$lgY[dat$non_i]
            F  <- Y %*% t(loglam_m) + cu - rep(q * Lm, each = nu)     # nu x M
            ls <- apply(F, 1, lse)
            nll <- nll - (sum(ls) - nu * log(dat$M))
            W  <- exp(F - ls)
            EL  <- (W %*% Lm)[, 1]
            EL2 <- (W %*% Lm^2)[, 1]
            SumEL <- SumEL + sum(EL)
            Vs    <- Vs + sum(EL2 - EL^2)
            Elam  <- Elam + (t(colSums(W)) %*% lam_m)[1, ]
        }

        ## collared animals with telemetry: Gauss-Hermite about the mean fix
        for (j in seq_along(dat$tel_i)) {
            i  <- dat$tel_i[j]; ti <- dat$telinfo[[j]]
            sdj <- sigma / sqrt(ti$nl)
            xg <- ti$tb[1] + sdj * dat$z1
            yg <- ti$tb[2] + sdj * dat$z2
            d2 <- (matrix(rep(xg, K), dat$G, K) - dat$Cm)^2 + (matrix(rep(yg, K), dat$G, K) - dat$Cy)^2
            loglam <- pars$log_lambda0 - d2 / two_s2
            lam    <- exp(loglam)
            Lg     <- (lam %*% U)[, 1]
            cu <- dat$YlogU[i] + dat$sumY[i] * log(q) - dat$lgY[i]
            F  <- (loglam %*% dat$CH[i, ])[, 1] + cu - q * Lg + dat$logw
            ls <- lse(F)
            tel <- -ti$nl * log(2 * pi * sigma^2) - ti$ssd / two_s2 + log(2 * pi * sigma^2 / ti$nl)
            nll <- nll - (tel - log(dat$A_mask) + ls)
            w   <- exp(F - ls)
            EL  <- sum(w * Lg); EL2 <- sum(w * Lg^2)
            SumEL <- SumEL + EL
            Vs    <- Vs + (EL2 - EL^2)
            Elam  <- Elam + (t(w) %*% lam)[1, ]
        }

        ## unidentified marked cues, by detector
        if (dat$has_Tm) {
            mean_k <- (1 - q) * U * Elam
            pos <- dat$Tm_pos
            nll <- nll - (sum(dat$Tm_k[pos] * log(mean_k[pos] + 1e-30)) - sum(mean_k) - dat$lgTm)
        }

        ## total density = unmarked + marked contribution
        D  <- Du + SumEL / J1

        ## overdispersion of the unmarked total, averaged over deployment
        mu <- Du * J1
        if (dat$distribution == "binomial") {
            r <- D * (J2 - J1^2 / dat$A_mask) / mu
        } else {
            r <- D * J2 / mu
        }
        chat <- if (dat$chat_on_tape) 1 + smooth_pos(r + Vs / mu, 1e-3) else dat$chat

        C  <- dat$C_u
        kd <- mu / (chat - 1)
        nll <- nll - (lgamma(C + kd) - lgamma(C + 1) - lgamma(kd) +
                      kd * log(kd / (kd + mu)) + C * log(mu / (kd + mu)))

        RTMB::REPORT(chat); RTMB::REPORT(mu); RTMB::REPORT(Vs); RTMB::REPORT(SumEL)
        RTMB::ADREPORT(D); RTMB::ADREPORT(lambda0); RTMB::ADREPORT(sigma)
        if (dat$has_Tm) RTMB::ADREPORT(q)
        nll
    }

    fit_once <- function (start) {
        obj <- RTMB::MakeADFun(nll_joint, start, map = map, silent = !trace)
        opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr,
                               control = list(iter.max = 500, eval.max = 1000)),
                        error = function(e) e)
        list(obj = obj, opt = opt)
    }

    ## 4. Two passes: c-hat on the tape, then held fixed ---------------------
    f1 <- fit_once(pars)
    if (inherits(f1$opt, "error")) return(NULL)
    chat_hat <- f1$obj$report(f1$obj$env$last.par.best)$chat
    dat$chat_on_tape <- FALSE
    dat$chat <- as.numeric(chat_hat)
    f <- fit_once(f1$obj$env$parList(f1$opt$par))
    if (inherits(f$opt, "error") || f$opt$convergence != 0) return(NULL)

    sdr <- RTMB::sdreport(f$obj)
    list(estimates = summary(sdr, select = "report"),
         chat      = dat$chat,
         has_Tm    = has_Tm,
         proctime  = unname((proc.time() - ptm)[3]))
}

###############################################################################
## extractfn for secrdesign::run.scenarios(fit = FALSE, ...)
##
## run.scenarios calls extractfn(capthist, ...) for each simulated capthist when
## fit = FALSE; the additional named arguments of run.scenarios are passed on, so supply
## the mask (and optionally 'distribution') there, e.g.
##   run.scenarios(nrepl, scenario, fit = FALSE, extractfn = rtmb_extract,
##                 mask = msk, distribution = "binomial", ...)
## Returns a data frame with one row per parameter (D, lambda0, sigma, plus pID if Tm are
## present) and columns estimate, SE.estimate, lcl, ucl, in the form returned for a
## secr fit by predict(), so that secrdesign::estimateSummary applies to the output.
## Intervals are log-normal as in secr. Attributes: proctime (seconds), chat.
## A failed fit gives NA estimates.
###############################################################################
rtmb_extract <- function (capthist, mask, distribution = "binomial", ...) {
    fit <- try(fit_rtmb_marking(capthist, mask, distribution), silent = TRUE)
    rows <- c("D", "lambda0", "sigma")
    if (inherits(fit, "try-error") || is.null(fit)) {
        out <- data.frame(estimate = NA_real_, SE.estimate = NA_real_, lcl = NA_real_, ucl = NA_real_,
                          row.names = "D")
        out <- out[rep(1, 3), ]
        rownames(out) <- rows
        attr(out, "proctime") <- NA_real_
        attr(out, "chat") <- NA_real_
        return(out)
    }
    est <- fit$estimates
    if (fit$has_Tm) rows <- c(rows, "q")
    est <- est[match(rows, rownames(est)), , drop = FALSE]
    z   <- qnorm(0.975)
    ## log-normal interval as secr: w = exp(z * sqrt(log(1 + (SE/estimate)^2)))
    w   <- exp(z * sqrt(log(1 + (est[, 2] / est[, 1])^2)))
    out <- data.frame(estimate = est[, 1], SE.estimate = est[, 2],
                      lcl = est[, 1] / w, ucl = est[, 1] * w, row.names = rows)
    if (fit$has_Tm) rownames(out)[rownames(out) == "q"] <- "pID"
    attr(out, "proctime") <- fit$proctime
    attr(out, "chat") <- fit$chat
    out
}

###############################################################################
## Standalone alternative to rtmb_extract for run.scenarios(fit = FALSE) with several masks.
## Instead of a single 'mask' it takes a list 'masks' and chooses the one for each scenario by
## the scenario's 'maskindex', found in the calling frame (the 'scenario' argument of
## secrdesign's internal processCH, which calls extractfn). All scenarios can then be run in
## one call:
##   run.scenarios(nrepl, scenarios, fit = FALSE, extractfn = rtmb_extract_masks,
##                 masks = list(mskW25, mskW50), maskset = list(mskW25, mskW50), ...)
## Otherwise the same as rtmb_extract.
###############################################################################
rtmb_extract_masks <- function (capthist, distribution = "binomial", ...) {
    mask <- get("fitarg", envir = parent.frame())$mask
    fit <- try(fit_rtmb_marking(capthist, mask, distribution), silent = TRUE)
    rows <- c("D", "lambda0", "sigma")
    if (inherits(fit, "try-error") || is.null(fit)) {
        out <- data.frame(estimate = NA_real_, SE.estimate = NA_real_, lcl = NA_real_, ucl = NA_real_,
                          row.names = "D")
        out <- out[rep(1, 3), ]
        rownames(out) <- rows
        attr(out, "proctime") <- NA_real_
        attr(out, "chat") <- NA_real_
        return(out)
    }
    est <- fit$estimates
    if (fit$has_Tm) rows <- c(rows, "q")
    est <- est[match(rows, rownames(est)), , drop = FALSE]
    z   <- qnorm(0.975)
    w   <- exp(z * sqrt(log(1 + (est[, 2] / est[, 1])^2)))    # log-normal interval as secr
    out <- data.frame(estimate = est[, 1], SE.estimate = est[, 2],
                      lcl = est[, 1] / w, ucl = est[, 1] * w, row.names = rows)
    if (fit$has_Tm) rownames(out)[rownames(out) == "q"] <- "pID"
    attr(out, "proctime") <- fit$proctime
    attr(out, "chat") <- fit$chat
    out
}

