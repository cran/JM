`opt.longWB` <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    w <- (logT - eta.t) / sigma.t
    ew <- - exp(w)
    log.p.tb <- d * (w - log(sigma.t)) + ew
    p.bytn <- p.byt * (log.p.yb + log.p.tb)
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}

