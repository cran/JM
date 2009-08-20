H.longWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.s <- alpha * Ys
    exp.eta.tw.P <- exp(eta.tw) * P
    H1 <- XtX / sigma^2
    Int <- wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s) * alpha^2
    H2 <- matrix(0, ncx, ncx)
    for (i in 1:ncx) {
        for (j in i:ncx) {
            ki <- exp.eta.tw.P * rowsum(Int * Xs[, i] * Xs[, j], id.GK, reorder = FALSE)
            kii <- c((p.byt * ki) %*% wGH)
            H2[i, j] <- sum(kii, na.rm = TRUE)
        }
    }
    H2[lower.tri(H2)] <- t(H2)[lower.tri(H2)]
    H1 + H2
}

