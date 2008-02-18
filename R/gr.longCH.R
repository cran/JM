`gr.longCH` <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    ew1 <- - exp(eta.t)
    out <- (d + c((ew1 * p.byt) %*% wGH)) * (alpha * Xtime)
    sc2 <- - colSums(out, na.rm = TRUE)
    c(sc1 + sc2)    
}

