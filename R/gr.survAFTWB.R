gr.survAFTWB <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    sigma.t <- if (is.null(scaleWB)) exp(thetas[ncww + 2]) else scaleWB
    eta.tw <- as.vector(WW %*% gammas)
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    wk.exp.eta.s <- wk * exp(eta.s)
    exp.eta.tw <- exp(eta.tw)
    Vi <- exp.eta.tw * P * rowsum(wk.exp.eta.s, id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    Vii <- d * (sigma.t - 1) / Vi - sigma.t * Vi^(sigma.t - 1)
    scgammas <- - colSums(WW * (d + c((p.byt * Vii * Vi) %*% wGH)), na.rm = TRUE)
    scalpha <- - sum((p.byt * (d * Y + Vii * exp.eta.tw * P * rowsum(wk.exp.eta.s * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    if (is.null(scaleWB)) {
        scsigmat <- - sigma.t * sum((p.byt * (d / sigma.t + (d - Vi^sigma.t) * log(Vi))) %*% wGH, na.rm = TRUE)
        c(scgammas, scalpha, scsigmat)
    } else {
        c(scgammas, scalpha)
    }
}

