opt.survAFTWB <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    sigma.t <- exp(thetas[ncww + 2])
    eta.tw <- as.vector(WW %*% gammas)
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
    log.survival <- - Vi^sigma.t
    log.p.tb <- d * log.hazard + log.survival    
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}

