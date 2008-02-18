`gr.longWB` <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.t <- c(WW %*% gammas) + Y * alpha
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    w <- (logT - eta.t) / sigma.t
    ki <- exp(w) - d
    kii <- c((p.byt * ki) %*% wGH)
    sc2 <- - colSums((alpha * Xtime) * kii, na.rm = TRUE) / sigma.t
    c(sc1 + sc2)
}

