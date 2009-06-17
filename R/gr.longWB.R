gr.longWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.s <- alpha * Ys
    exp.eta.tw.P <- exp(eta.tw) * P
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Int <- wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s) * alpha
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- exp.eta.tw.P * rowsum(Int * Xs[, i], id.GK, reorder = FALSE)
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- - sum(d * alpha * Xtime[, i] - kii, na.rm = TRUE)
    }    
    c(sc1 + sc2)
}

