`update.bCH` <-
function (b, hes.b, betas, sigma, thetas, D, transformed = FALSE) {
    log.p.yt <- numeric(n)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    gammas <- thetas[-(ncww + 1)]
    alpha <- thetas[ncww + 1]
    if (transformed)
        gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
    eta.tw <- as.vector(WW %*% gammas)
    sc <- as.vector(S %*% diff(gammas[1:nk]))
    environment(fn.b) <- environment(gr.b) <- environment(logsurvCH) <- environment()
    for (i in 1:n) {
        # individual values
        ind.i <- id == i
        yi <- y[ind.i]
        eta.yxi <- eta.yx[ind.i]
        yi.eta.yxi <- yi - eta.yxi
        logTi <- logT[i]
        eta.yxTi <- eta.yxT[i]
        eta.twi <- eta.tw[i]
        Z.ind.i <- Z[ind.i, , drop = FALSE]
        Ztime.i <- Ztime[i, ]
        # posterior modes
        opt <- try(optim(b[i, ], fn.b, gr.b, method = "BFGS", hessian = TRUE), TRUE)
        opt <- if (!inherits(opt, "try-error")) {
            opt 
        } else {
           list(par = b[i, ], hessian = matrix(hes.b[i, ], ncz, ncz))
        }
        H <- opt$hessian
        var.b <- if (!inherits(var.b <- try(solve(H), TRUE), "try-error")) var.b else ginv(H)
        b[i, ] <- opt$par
        vb[i, ] <- c(var.b)
        hes.b[i, ] <- c(H)
        # likelihood contributions
        mu.y.b <- eta.yxi + rowSums(Z.ind.i * rep(b[i, ], each = ni[i]))
        log.p.y.b <- sum(dnorm(yi, mu.y.b, sigma, log = TRUE))
        log.p.t.b <- logsurvCH(alpha, b[i, ])
        log.p.b <- if (diag.D) {
            sum(dnorm(b[i, ], 0, sqrt(D), log = TRUE))
        } else {
            dmvnorm(b[i, ], rep(0, ncz), D, log = TRUE)
        }
        log.p.yt[i] <- (log.p.y.b + log.p.t.b + log.p.b) - 0.5 * log(det(H))
    }
    res <- cons.logLik + sum(log.p.yt, na.rm = TRUE)
    attr(res, "b") <- b
    attr(res, "vb") <- vb
    attr(res, "hes.b") <- hes.b
    res
}

