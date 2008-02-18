`opt.survWB` <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    sigma.t <- exp(thetas[ncww + 2])
    eta <- c(WW %*% gammas) + Y * alpha
    w <- (logT - eta) / sigma.t
    ew <- - exp(w)
    log.p.tb <- d * (w - log(sigma.t)) + ew
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}

