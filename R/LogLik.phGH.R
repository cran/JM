`LogLik.phGH` <-
function (thetas) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    alpha <- thetas[ncx + ncww + 2]
    D <- thetas[seq(ncx + ncww + 3, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.yxT2 <- as.vector(Xtime2 %*% betas)
    Y <- eta.yxT + Ztime.b
    Y2 <- eta.yxT2 + Ztime2.b
    eta.t <- if (is.null(WW)) alpha * Y else as.vector(WW %*% gammas) + alpha * Y
    eta.t2 <- if (is.null(WW)) alpha * Y2 else as.vector(WW %*% gammas)[indT] + alpha * Y2    
    lambda0T <- d * lambda0[ind.T0]
    lambda0T[is.na(lambda0T)] <- 0
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    ew <- exp(eta.t2)
    log.p.tb <- d * (log(lambda0T) + eta.t)
    log.p.tb[ind0, ] <- 0
    log.p.tb[ind.lenN0, ] <- log.p.tb[ind.lenN0, ] - rowsum(lambda0[ind.lambda] * ew, indT)
    log.p.b <- if (ncz == 1) {
        dnorm(b, sd = sqrt(D), log = TRUE)
    } else {
        if (diag.D) {
            rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
        } else {
            dmvnorm(b, rep(0, ncz), D, TRUE)
        }
    }
    p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n))
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    p.byt <- p.ytb / p.yt
    log.p.yt <- log(p.yt)
    - sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
}

