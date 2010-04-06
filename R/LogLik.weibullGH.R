LogLik.weibullGH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    sigma.t <- if (is.null(scaleWB)) exp(thetas$log.sigma.t) else scaleWB
    D <- thetas$D
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    Y <- eta.yxT + Ztime.b
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)    
    log.hazard <- log(sigma.t) + (sigma.t - 1) * logT + eta.t
    log.survival <- - exp(eta.tw) * P * rowsum(wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- d * log.hazard + log.survival
    log.p.b <- if (ncz == 1) {
        dnorm(b, sd = sqrt(D), log = TRUE)
    } else {
        if (diag.D) {
            rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
        } else {
            dmvnorm(b, rep(0, ncz), D, TRUE)
        }
    }
    p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n)); dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    log.p.yt <- log(p.yt)
    - sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
}

