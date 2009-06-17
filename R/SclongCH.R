SclongCH <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    mu <- y - eta.yx
    sc1 <- - crossprod(X, mu - Zb) / sigma^2
    Y <- eta.yxT + rowSums(Ztime * b)
    eta.t <- eta.tw + Y * alpha
    ew1 <- - exp(eta.t)
    ew2 <- ew1 - trc.t1
    out <- (d + ew2) * (alpha * Xtime)
    sc2 <- - colSums(out, na.rm = TRUE)
    c(sc1 + sc2)    
}

