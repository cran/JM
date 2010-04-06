gr.survWB <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    sigma.t <- if (is.null(scaleWB)) exp(thetas[ncww + 2]) else scaleWB
    eta.tw <- as.vector(WW %*% gammas)
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    exp.eta.tw.P <- exp(eta.tw) * P
    Int <- wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s)
    ki <- exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
    kii <- c((p.byt * ki) %*% wGH)
    scgammas <- - colSums(WW * (d - kii), na.rm = TRUE)
    scalpha <- - sum((p.byt * (d * Y - exp.eta.tw.P * rowsum(Int * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)    
    Int2 <- st^(sigma.t - 1) * (1 + sigma.t * log.st) * exp(eta.s)
    if (is.null(scaleWB)) {
        scsigmat <- - sigma.t * sum((p.byt * (d * (1/sigma.t + logT) - exp.eta.tw.P * rowsum(wk * Int2, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
        c(scgammas, scalpha, scsigmat)
    } else {
        c(scgammas, scalpha)
    }
}

