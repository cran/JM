gr.longPC <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.s <- alpha * Ys
    exp.eta.tw <- exp(eta.tw)
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Int <- xi[ind.K] * wkP * exp(eta.s) * alpha
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- exp.eta.tw * rowsum(Int * Xs[, i], id.GK, reorder = FALSE)
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- - sum(d * alpha * Xtime[, i] - kii, na.rm = TRUE)
    }    
    c(sc1 + sc2)
}

