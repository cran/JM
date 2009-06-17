gr.longAFTWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.s <- alpha * Ys
    wk.exp.eta.s <- wk * exp(eta.s)
    exp.eta.tw.P <- exp(eta.tw) * P
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Vi <- exp(eta.tw) * P * rowsum(wk.exp.eta.s, id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    Vii <- d * (sigma.t - 1) / Vi - sigma.t * Vi^(sigma.t - 1)
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- Vii * exp.eta.tw.P * rowsum(wk.exp.eta.s * alpha * Xs[, i], id.GK, reorder = FALSE)
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- - sum(d * alpha * Xtime[, i] + kii, na.rm = TRUE)
    }    
    c(sc1 + sc2)
}

