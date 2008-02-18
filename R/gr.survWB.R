`gr.survWB` <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    sigma.t <- exp(thetas[ncww + 2])
    eta <- c(WW %*% gammas) + Y * alpha
    w <- (logT - eta) / sigma.t
    ki <- exp(w) - d
    kmi <- w * ki
    kii <- c((p.byt * ki) %*% wGH)
    scgammas <- - colSums(WW * kii, na.rm = TRUE) / sigma.t
    scalpha <- - sum((p.byt * Y * ki) %*% wGH, na.rm = TRUE) / sigma.t
    scsigmat <- - sigma.t * sum(c((p.byt * kmi) %*% wGH) - d, na.rm = TRUE) / sigma.t
    c(scgammas, scalpha, scsigmat)
}

