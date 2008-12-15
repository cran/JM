`posterior.b` <-
function (b) {
    Y <- eta.yxT + rowSums(Ztime.missO * b)
    eta.t <- eta.tw + alpha.new * Y
    mu.y <- eta.yx + rowSums(Z.missO * b[id3.miss, , drop = FALSE])
    logNorm <- dnorm(y.missO, mu.y, sigma.new, TRUE)
    log.p.yb <- as.vector(tapply(logNorm, id.miss, sum))
    log.p.tb <- if (method == "weibull-GH") {
        w <- (logT.missO - eta.t) / sigma.t.new
        ew <- - exp(w)
        d.missO * (w - log(sigma.t.new)) + ew
    } else {
        kn <- object$knots
        ord <- object$control$ord
        S <- splineDesign(kn[-c(1, length(kn))], logT.missO, ord = ord - 1)
        S <- ord * S / rep(diff(kn, lag = ord + 1), each = length(d.missO))
        sc <- as.vector(S %*% diff(gammas.new[1:nk]))
        ew <- - exp(eta.t)
        d.missO * (log(sc) + eta.t - logT.missO) + ew
    }
    log.p.b <- if (ncz == 1) {
        dnorm(b, sd = sqrt(D.new), log = TRUE)
    } else {
        if (diag.D) {
            rowSums(dnorm(b, sd = rep(sqrt(D.new), each = k), log = TRUE))
        } else {
            dmvnorm(b, rep(0, ncz), D, TRUE)
        }
    }
    log.p.yb + log.p.tb + log.p.b
}

