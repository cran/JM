LogLik.piecewiseGH <-
function (thetas) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    alpha <- thetas[ncx + ncww + 2]
    xi <- exp(thetas[seq(ncx + ncww + 3, ncx + ncww + 2 + Q)])
    D <- thetas[seq(ncx + ncww + Q + 3, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)    
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
    Y <- eta.yxT + Ztime.b
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)    
    log.hazard <- log(xi[ind.D]) + eta.t
    log.survival <- - exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
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

