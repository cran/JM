`posterior.b` <-
function (b) {
    Y <- eta.yxT + rowSums(Ztime.missO * b)
    eta.t <- eta.tw + alpha.new * Y
    mu.y <- eta.yx + rowSums(Z.missO * b[id3.miss, , drop = FALSE])
    logNorm <- dnorm(y.missO, mu.y, sigma.new, TRUE)
    log.p.yb <- as.vector(tapply(logNorm, id.miss, sum))
    log.p.tb <- if (method == "ph-GH") {
        indT <- object$indexes$indT # change to indT.O
        ind.lambda <- object$indexes$ind.lambda # change to ind.lambda.O
        lambda0 <- object$coefficients$lambda0[, "basehaz"]
        Time <- exp(logT) # change to Time.O
        unqT <- sort(unique(Time[d == 1])) # change to unqT.O
        ind.T0 <- match(Time, unqT) # change to lambda0T.O
        ind.T0[d == 0] <- NA # change to lambda0T.O
        ind.lenN0 <- object$indexes$ind.len != 0 # change to ind.lenN0.O
        lambda0T <- lambda0[ind.T0] # change to lambda0T.O
        lambda0T[is.na(lambda0T)] <- 0
        eta.yxT2 <- as.vector(Xtime2 %*% betas.new) # change Xtime2 -> Xtime2.missO
        Y2 <- eta.yxT2 + rowSums(Ztime2 * b[id4.miss, ]) # change Ztime2 -> Ztime2.missO; compute id4.miss
        eta.t2 <- if (is.null(WW)) alpha.new * Y2 else as.vector(WW.missO %*% gammas.new)[indT] + alpha.new * Y2
        ew <- exp(eta.t2)
        out <- d.missO * (log(lambda0T) + eta.t)
        out[ind0, ] <- 0
        out[ind.lenN0, ] <- out[ind.lenN0, ] - rowsum(lambda0[ind.lambda] * ew, indT)
        out
    } else if (method == "weibull-GH") {
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

