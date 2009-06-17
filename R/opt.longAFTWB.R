opt.longAFTWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
    log.survival <- - Vi^sigma.t
    log.p.tb <- d * log.hazard + log.survival
    p.bytn <- p.byt * (log.p.yb + log.p.tb)
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}

