opt.survWB <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    sigma.t <- exp(thetas[ncww + 2])
    eta.tw <- as.vector(WW %*% gammas)
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    log.hazard <- log(sigma.t) + (sigma.t - 1) * logT + eta.t
    log.survival <- - exp(eta.tw) * P * rowsum(wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- d * log.hazard + log.survival    
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}

