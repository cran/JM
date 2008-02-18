`LogLik.weibullGH` <-
function (thetas) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    alpha <- thetas[ncx + ncww + 2]
    sigma.t <- exp(thetas[ncx + ncww + 3])
    D <- thetas[seq(ncx + ncww + 4, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    w <- (logT - eta.t) / sigma.t
    ew <- - exp(w)
    log.p.tb <- d * (w - log(sigma.t)) + ew
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

