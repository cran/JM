gr.survAFTWB <-
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
    wk.exp.eta.s <- wk * exp(eta.s)
    exp.eta.tw <- exp(eta.tw)
    Vi <- exp.eta.tw * P * rowsum(wk.exp.eta.s, id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    Vii <- d * (sigma.t - 1) / Vi - sigma.t * Vi^(sigma.t - 1)
    scgammas <- - colSums(WW * (d + c((p.byt * Vii * Vi) %*% wGH)), na.rm = TRUE)
    scalpha <- if (parameterization %in% c("value", "both")) {
        - sum((p.byt * (d * Y + Vii * exp.eta.tw * P * rowsum(wk.exp.eta.s * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        - sum((p.byt * (d * Y.deriv + Vii * exp.eta.tw * P * rowsum(wk.exp.eta.s * Ys.deriv, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    scsigmat <- if (is.null(scaleWB)) {
         - sigma.t * sum((p.byt * (d / sigma.t + (d - Vi^sigma.t) * log(Vi))) %*% wGH, na.rm = TRUE)
        
    } else NULL
    c(scgammas, scalpha, scalpha.D, scsigmat)
}

