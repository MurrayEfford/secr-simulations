# LDF material from secr-simulation.rmd
# put aside 2024-03-25
# - scenario of Ben Stevenson retained

```{r LDF old shared, eval = TRUE}
detarg <- expand.detarg(seq(0,1.6,0.2), 150, cov.structure = "sq_exponential",
                        noccasions = 1, shared = TRUE)
poparg <- list(buffer = 500)
fitarg <- list(detectfn = 'HHN', mask = msk)
scen   <- make.scenarios(D = 1, lambda0 = 0.5, sigma = 100,  detectfn = 'HHN', 
                         detindex = 1:length(detarg))
```

```{r LDF old shared run, eval = FALSE}
LDFsimsoldshared <- run.scenarios(500, scen, trapset = trps, CH.function = "simLDF.capthist",
                                  pop.args = poparg, det.args = detarg, fit = TRUE, fit.args = fitarg)
saveRDS(LDFsimsoldshared, file = 'LDFsimsoldshared.RDS')
```

```{r LDFoldfigs, echo = FALSE}
LDFsimsold <- readRDS(file = 'LDFsimsold.RDS')
LDFsimsoldshared <- readRDS(file = 'LDFsimsoldshared.RDS')
estDold    <- estimateSummary(LDFsimsold,'D')
estDoldshared <- estimateSummary(LDFsimsoldshared,'D', validrange = c(0,10))
estLold    <- estimateSummary(LDFsimsold,'lambda0')
estLoldshared <- estimateSummary(LDFsimsoldshared,'lambda0', validrange = c(0,5))
estSold    <- estimateSummary(LDFsimsold,'sigma')
estSoldshared <- estimateSummary(LDFsimsoldshared,'sigma', validrange = c(0,500))
sigma.u <- sapply(detarg, function(x) x$cov.pars$sigma.u)
cvu <- sqrt(exp(sigma.u^2)-1)
leg <- c('D','lambda0','sigma')
par(mfrow=c(1,2))

plot(0,0, type = 'n', xlim=c(0,1.6), xlab = 'sigma.u', ylab = 'RB', ylim=c(-0.4,0.8))
addRB(sigma.u, estLoldshared, pch=21, bg='white')
addRB(sigma.u, estSoldshared, pch=24, bg='white')
addRB(sigma.u, estDoldshared, pch=16)
abline(h=0, lty=2)
legend(0.05,0.8, legend=leg, pch=c(16,21,24))
mtext(side=3, 'Common LDF', cex=0.9)

plot(0,0, type = 'n', xlim=c(0,1.6), xlab = 'sigma.u', ylab = 'RB', ylim=c(-0.4,0.8))
addRB(sigma.u, estLold, pch=21, bg='white')
addRB(sigma.u, estSold, pch=24, bg='white')
addRB(sigma.u, estDold, pch=16)
abline(h=0, lty=2)
legend(0.05,0.8, legend=leg, pch=c(16,21,24))
mtext(side=3, 'Individual LDF', cex=0.9)
```

```{r LDF, eval =FALSE}
library (secrdesign)
nc <- setNumThreads(18)
grid <- make.grid(nx = 8, ny = 8, spacing = 40, detector = "count")
detarg <- list(cov.structure = "sq_exponential", 
               cov.pars = list(sigma.u = 0.5, rho = 40))
detarg <- rep(list(detarg), 10)
#for (i in 1:10) detarg[[i]]$cov.pars$rho <- i * 20
for (i in 1:10) detarg[[i]]$cov.pars$sigma.u <- (i-1) * 0.2

fitarg <- list(detectfn = 'HHN')
scen <- make.scenarios(D=10, lambda0 = 0.5, sigma = 20,  detectfn = 'HHN', detindex = 1:10)
LDFsims <- run.scenarios(500, scen, trapset = grid, CH.function = "simLDF.capthist",
                         det.args = detarg, fit = TRUE, fit.args = fitarg)
saveRDS(LDFsims, file = 'LDFsims.RDS')
```

## Results

```{r LDF sumamry}
LDFsims <- readRDS(file = 'LDFsims.RDS')
counts <- t(sapply(LDFsims$output, function(x) apply(sapply(x, attr, 'counts'),1,mean)))
sigma.u <- seq(0,1.8,0.2)
estD <- estimateSummary(LDFsims)
estL <- estimateSummary(LDFsims, 'lambda0')
estS <- estimateSummary(LDFsims, 'sigma')

par(pty='s')
plot(0,0,type='n', xlim=c(0,2), ylim=c(-0.4,0.4), xlab = 'sigma.u', ylab = 'RB')
addRB(sigma.u, estD, pch=16, type='o')
addRB(sigma.u, estL, pch=21, bg = 'white', type='o')
addRB(sigma.u, estS, pch=24, bg = 'white', type='o')
abline (h=0, lty=2)
```
