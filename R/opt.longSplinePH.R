opt.longSplinePH <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.t <- eta.tw2 + eta.tw1 + alpha * Y
    eta.s <- alpha * Ys
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    log.hazard <- eta.t
    log.survival <- - exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE)
    log.p.tb <- d * log.hazard + log.survival
    p.bytn <- p.byt * (log.p.yb + log.p.tb)
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}

