`fn.b` <-
function (bb) {
    eta.yi <- eta.yxi + rowSums(Z.ind.i * rep(bb, each = ni[i]))
    log.p.ybi <- sum(dnorm(yi, eta.yi, sigma, log = TRUE))
    eta.ti <- eta.twi + alpha * (eta.yxT[i] + sum(Ztime.i * bb))
    log.p.tbi <- if (d[i]) log(sc[i]) + eta.ti - exp(eta.ti) - logT[i] else - exp(eta.ti)
    log.p.bi <- if (diag.D) - 0.5 * crossprod(bb, bb / D)[1, ] else - 0.5 * crossprod(bb, solve(D, bb))[1, ]
    - (log.p.ybi + log.p.tbi + log.p.bi)
}

