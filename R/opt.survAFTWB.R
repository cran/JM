opt.survAFTWB <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    sigma.t <- if (is.null(scaleWB)) exp(thetas$log.sigma.t) else scaleWB
    eta.tw <- as.vector(WW %*% gammas) 
    eta.t <- switch(parameterization, "value" = eta.tw + alpha * Y, 
        "slope" = eta.tw + Dalpha * Y.deriv, "both" = eta.tw + alpha * Y + Dalpha * Y.deriv)    
    eta.s <- switch(parameterization, "value" = alpha * Ys, "slope" = Dalpha * Ys.deriv, 
        "both" = alpha * Ys + Dalpha * Ys.deriv)
    Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
    log.survival <- - Vi^sigma.t
    log.p.tb <- d * log.hazard + log.survival    
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}

