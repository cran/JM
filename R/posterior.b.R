posterior.b <-
function (b) {
    Y <- eta.yxT + rowSums(Ztime.missO * b)
    eta.t <- eta.tw + alpha.new * Y
    mu.y <- eta.yx + rowSums(Z.missO * b[id3.miss, , drop = FALSE])
    logNorm <- dnorm(y.missO, mu.y, sigma.new, TRUE)
    log.p.yb <- as.vector(tapply(logNorm, id.miss, sum))
    log.p.tb <- if (method == "weibull-PH-GH") {
        id.GK3 <- rep(seq_len(n.missO), each = object$control$GKk)
        Ys.new <- as.vector(Xs.missO %*% betas.new) +rowSums(Zs.missO * b[id.GK3, , drop = FALSE])
        eta.s <- alpha.new * Ys.new
        wk <- rep(wk, n.missO)
        log.hazard <- log(sigma.t.new) + (sigma.t.new - 1) * logT.missO + eta.t
        Vi <- exp(log(sigma.t.new) + (sigma.t.new - 1) * log.st.missO + eta.s)
        log.survival <- - exp(eta.tw) * P.missO * as.vector(tapply(wk * Vi, id.GK3, sum))
        d.missO * log.hazard + log.survival
    } else if (method == "weibull-AFT-GH") {
        id.GK3 <- rep(seq_len(n.missO), each = object$control$GKk)
        Ys.new <- as.vector(Xs.missO %*% betas.new) +rowSums(Zs.missO * b[id.GK3, , drop = FALSE])
        eta.s <- alpha.new * Ys.new
        wk <- rep(wk, n.missO)
        Vi <- exp(eta.tw) * P.missO * as.vector(tapply(wk * exp(eta.s), id.GK3, sum))
        log.hazard <- log(sigma.t.new) + (sigma.t.new - 1) * log(Vi) + eta.t
        log.survival <- - Vi^sigma.t.new
        d.missO * log.hazard + log.survival
    } else if (method == "spline-PH-GH") {
        id.GK3 <- rep(seq_len(n.missO), each = object$control$GKk)
        Ys.new <- as.vector(Xs.missO %*% betas.new) +rowSums(Zs.missO * b[id.GK3, , drop = FALSE])
        eta.s <- alpha.new * Ys.new
        wk <- rep(wk, n.missO)
        log.hazard <- c(W2.missO %*% gammas.bs.new) + eta.t
        Vi <- exp(c(W2s.missO %*% gammas.bs.new) + eta.s)
        log.survival <- - exp(eta.tw) * P.missO * as.vector(tapply(wk * Vi, id.GK3, sum))
        d.missO * log.hazard + log.survival
    } else if (method == "piecewise-PH-GH") {
        ii <- object$x$id.GK[id.GK]
        nn <- as.vector(tapply(ii, ii, length))
        id.GK3 <- rep(seq_len(n.missO), nn)
        Ys.new <- as.vector(Xs.missO %*% betas.new) +rowSums(Zs.missO * b[id.GK3, , drop = FALSE])
        eta.s <- alpha.new * Ys.new
        log.hazard <- log(xi.new[ind.D.missO]) + eta.t
        log.survival <- - exp(eta.tw) * as.vector(tapply(xi.new[ind.K.missO] * wkP.missO * exp(eta.s), id.GK3, sum))
        d.missO * log.hazard + log.survival
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

