`gr.longPH` <-
function (X, Xtime, Xtime2, ew) {
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    nw <- ncol(Xtime2)
    res <- numeric(nw)
    for (i in 1:nw) {
        p <- rowsum(lambda0. * (alpha * Xtime2[, i]) * ew, indT)
        res[i] <- sum((p * p.byt.) %*% wGH)
    }
    sc2 <- - (crossprod(alpha * Xtime, d) - res)
    c(sc1 + sc2)
}

